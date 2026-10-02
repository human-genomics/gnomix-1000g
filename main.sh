#!/usr/bin/env bash
# =============================================================================
# gnomix-1000g — local ancestry for the 3,202 high-coverage 1000 Genomes samples
# =============================================================================
# Single entry point. Docker needs nothing else. Local runs need bash and curl
# on Linux x86_64; the script installs Python 3.9 and Java if they are missing.
#
# Usage:
#   bash main.sh            # all steps
#   bash main.sh tracts     # one step (prepare, infer, trio, tracts, karyograms,
#                           #   compare, pfile, report; release packages outputs)
#
# Settings (environment variables, see README):
#   THREADS=4            worker processes (up to ~2.5 GB RAM each)
#   BATCH_SIZE=16        samples per inference task
#   CHROMS="1-22"        autosomes to process
#   TRIO_CHROMS="18-22"  chromosomes for the trio phasing benchmark
#   TRIO_GROUPS=6        Beagle runs per benchmark chromosome
#   BEAGLE_MEM_GB=16     Java heap for Beagle
#   RUN_TRIO_BENCHMARK=0 1: repeat the trio benchmark that decides where Gnofix is
#                        used (slowest step). 0: use the published decision in
#                        reference/gnofix_policy.json (same data, same result)
#   MAKE_PFILE=1         write the phase-corrected panel (output/pfile, ~10 GB)
#   GNOMIX1000G_DATA_DIR where downloads/ work/ output/ logs/ go (default: this directory)
#   REDRAW=1             redraw all karyograms
#
# Steps:
#   1. Tools: Python 3.9 environment, Gnomix (pinned commit), Beagle 5.5
#   2. Download inputs (6.4 GB; SHA-256 checked)
#   3. prepare     unpack models, index the panel, lift and match model SNPs
#   4. infer       pretrained Gnomix with and without Gnofix, every sample
#   5. trio        Gnofix policy per ancestry group (published, or re-run the benchmark)
#   6. tracts      msp files, haplotype tracts, global ancestry
#   7. karyograms  one PNG per sample
#   8. compare     Martin et al. (2017) and 1000 Genomes Phase 1 LAI
#   9. pfile       panel with the phase of the final tracts (Gnofix where used)
#  10. report      output/REPORT.md
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${GNOMIX1000G_DATA_DIR:-${SCRIPT_DIR}}"
STEP="${1:-all}"
export GNOMIX1000G_DATA_DIR="${DATA_DIR}"
export PYTHONPATH="${SCRIPT_DIR}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONWARNINGS=ignore PYTHONUNBUFFERED=1
export MPLCONFIGDIR="${MPLCONFIGDIR:-${DATA_DIR}/work/.matplotlib}"

# shellcheck source=lib.sh
source "${SCRIPT_DIR}/lib.sh"

# ---- 1. Tools ------------------------------------------------------------------
# The Docker image sets GNOMIX1000G_PYTHON and has its tools in /opt/tools.
PYTHON="${GNOMIX1000G_PYTHON:-}"
if [[ -z "${PYTHON}" ]]; then
    # Local run: tools/ inside the repository.
    export GNOMIX1000G_TOOLS="${GNOMIX1000G_TOOLS:-${SCRIPT_DIR}/tools}"
    TOOLS="${GNOMIX1000G_TOOLS}"
    mkdir -p "${TOOLS}"
    VENV="${TOOLS}/venv"
    if [[ ! -x "${VENV}/bin/python3" ]]; then
        BASE_PY=""
        for cand in python3.9 python3.10; do
            if command -v "${cand}" >/dev/null 2>&1; then BASE_PY="$(command -v "${cand}")"; break; fi
        done
        if [[ -z "${BASE_PY}" ]]; then
            if [[ "$(uname -s)-$(uname -m)" != "Linux-x86_64" ]]; then
                echo "ERROR: Python 3.9 or 3.10 is required (the pretrained models need" >&2
                echo "  scikit-learn 1.0.1 and xgboost 1.1.1). Install one, or use Docker." >&2
                exit 1
            fi
            echo "[setup] Installing standalone Python 3.9 into tools/python ..."
            fetch "https://github.com/astral-sh/python-build-standalone/releases/download/20241016/cpython-3.9.20%2B20241016-x86_64-unknown-linux-gnu-install_only.tar.gz" \
                "${TOOLS}/python.tar.gz" c20ee831f7f46c58fa57919b75a40eb2b6a31e03fd29aaa4e8dab4b9c4b60d5d
            tar -xzf "${TOOLS}/python.tar.gz" -C "${TOOLS}"
            BASE_PY="${TOOLS}/python/bin/python3.9"
        fi
        echo "[setup] Creating Python venv in tools/venv with ${BASE_PY} ..."
        "${BASE_PY}" -m venv "${VENV}"
    fi
    # (Re)install the pinned packages unless this exact requirements.txt installed successfully.
    REQ_HASH="$(sha256 "${SCRIPT_DIR}/requirements.txt")"
    if [[ "$(cat "${VENV}/.requirements.sha256" 2>/dev/null)" != "${REQ_HASH}" ]]; then
        echo "[setup] Installing pinned Python packages ..."
        "${VENV}/bin/python3" -m pip install --quiet --upgrade "pip<25"
        "${VENV}/bin/python3" -m pip install --quiet --no-cache-dir -r "${SCRIPT_DIR}/requirements.txt"
        echo "${REQ_HASH}" > "${VENV}/.requirements.sha256"
    fi
    PYTHON="${VENV}/bin/python3"
    # xgboost needs the OpenMP runtime (libgomp.so.1). Most systems have it;
    # if not, use a pinned copy from conda-forge (GCC runtime library exception).
    if ! "${PYTHON}" -c "import ctypes.util, sys; sys.exit(0 if ctypes.util.find_library('gomp') else 1)"; then
        if [[ ! -f "${TOOLS}/lib/libgomp.so.1" ]]; then
            echo "[setup] Installing OpenMP runtime (libgomp) into tools/lib ..."
            fetch "https://conda.anaconda.org/conda-forge/linux-64/libgomp-12.2.0-h65d4601_19.tar.bz2" \
                "${TOOLS}/libgomp.tar.bz2" 81a76d20cfdee9fe0728b93ef057ba93494fd1450d42bc3717af4e468235661e
            "${PYTHON}" -c "import tarfile, sys; tarfile.open(sys.argv[1]).extractall(sys.argv[2], members=[m for m in tarfile.open(sys.argv[1]) if m.name.startswith('lib/')])" \
                "${TOOLS}/libgomp.tar.bz2" "${TOOLS}"
            ln -sf libgomp.so.1.0.0 "${TOOLS}/lib/libgomp.so.1"
        fi
        export LD_LIBRARY_PATH="${TOOLS}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    fi
    if ! command -v java >/dev/null 2>&1; then
        if [[ ! -x "${TOOLS}/jre/bin/java" ]]; then
            echo "[setup] Installing Java 17 runtime into tools/jre ..."
            fetch "https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.20.1%2B1/OpenJDK17U-jre_x64_linux_hotspot_17.0.20.1_1.tar.gz" \
                "${TOOLS}/jre.tar.gz" 0b2b640e3046b64c8ec504de0ab9d91bb5610182bda21fad454681ce54d45a62
            mkdir -p "${TOOLS}/jre"
            tar -xzf "${TOOLS}/jre.tar.gz" -C "${TOOLS}/jre" --strip-components=1
        fi
        export PATH="${TOOLS}/jre/bin:${PATH}"
    fi
    bash "${SCRIPT_DIR}/setup_tools.sh" "${TOOLS}"
fi

# ---- Logging -------------------------------------------------------------------
if ! mkdir -p "${DATA_DIR}/logs" "${DATA_DIR}/output" 2>/dev/null; then
    echo "ERROR: cannot write to ${DATA_DIR}." >&2
    echo "  Docker with --user: create the host folder before 'docker run' (mkdir -p gnomix-1000g-data)." >&2
    exit 1
fi
mkdir -p "${MPLCONFIGDIR}"
LOG_FILE="${DATA_DIR}/logs/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1
echo "gnomix-1000g: step=${STEP} DATA_DIR=${DATA_DIR}"
echo "  THREADS=${THREADS:-4} CHROMS=${CHROMS:-1-22} RUN_TRIO_BENCHMARK=${RUN_TRIO_BENCHMARK:-0} TRIO_CHROMS=${TRIO_CHROMS:-18-22}"
echo "  MAKE_PFILE=${MAKE_PFILE:-1} REDRAW=${REDRAW:-0} BATCH_SIZE=${BATCH_SIZE:-16} TRIO_GROUPS=${TRIO_GROUPS:-6} BEAGLE_MEM_GB=${BEAGLE_MEM_GB:-16}"

# ---- 2. Downloads ----------------------------------------------------------------
if [[ "${STEP}" == "all" || "${STEP}" == "prepare" ]]; then
    echo "==> download inputs"
    bash "${SCRIPT_DIR}/download.sh" "${DATA_DIR}/downloads"
fi

# ---- 3-10. Pipeline ----------------------------------------------------------------
"${PYTHON}" -m gnomix1000g "${STEP}"
echo "Done. Outputs are in ${DATA_DIR}/output (start with output/REPORT.md)."
