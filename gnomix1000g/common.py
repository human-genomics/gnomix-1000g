"""Paths, constants and small helpers shared by all pipeline steps."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from contextlib import contextmanager
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("GNOMIX1000G_DATA_DIR", REPO_DIR))
TOOLS_DIR = Path(os.environ.get("GNOMIX1000G_TOOLS", REPO_DIR / "tools"))

DOWNLOADS = DATA_DIR / "downloads"
WORK = DATA_DIR / "work"
OUTPUT = DATA_DIR / "output"
REFERENCE = REPO_DIR / "reference"

AUTOSOMES = list(range(1, 23))

# Gnomix pretrained model labels, in model order (checked against each model).
ANCESTRIES = ["EUR", "EAS", "NAT", "AFR", "SAS", "AHG", "OCE", "WAS"]
# Names as in the Gnomix paper (Hilmarsson et al., Nat Commun 2026).
ANCESTRY_NAMES = {
    "EUR": "European",
    "EAS": "East Asian",
    "NAT": "Indigenous American",
    "AFR": "African",
    "SAS": "South Asian",
    "AHG": "African Hunter-Gatherer",
    "OCE": "Oceanian",
    "WAS": "West Asian",
}
# Colour per ancestry: hues of a validated categorical palette, assigned so that
# ancestries that occur together in one genome stay apart (checked with a CVD /
# normal-vision palette validator for each population context).
ANCESTRY_COLORS = {
    "AFR": "#eda100",
    "AHG": "#e87ba4",
    "EUR": "#2a78d6",
    "WAS": "#1baf7a",
    "SAS": "#008300",
    "EAS": "#4a3aa7",
    "NAT": "#e34948",
    "OCE": "#eb6834",
    "UNK": "#BDBDBD",
}

# Admixed populations. Their karyograms show EUR and WAS as one "European + West Asian"
# category, and the trio benchmark groups them apart from their super-population.
ADMIXED_POPULATIONS = ["ACB", "ASW", "CLM", "MXL", "PEL", "PUR"]
GROUPS = {"ACB": "AFR-American", "ASW": "AFR-American",
          "CLM": "AMR", "MXL": "AMR", "PEL": "AMR", "PUR": "AMR"}
# Karyogram colours of the merged categories of admixed populations.
WEST_EURASIAN = ANCESTRY_COLORS["EUR"]
OTHER = ANCESTRY_COLORS["EAS"]
# Continental labels: the closely related ancestries merged (AFR + AHG, EUR + WAS).
CONTINENTAL = {"AFR": "AFR", "AHG": "AFR", "EUR": "EUR", "WAS": "EUR",
               "NAT": "NAT", "EAS": "EAS", "SAS": "SAS", "OCE": "OCE"}

# 1000 Genomes population descriptions (IGSR).
POPULATION_NAMES = {
    "ACB": "African Caribbean in Barbados", "ASW": "African Ancestry in SW USA",
    "BEB": "Bengali in Bangladesh", "CDX": "Chinese Dai in Xishuangbanna, China",
    "CEU": "Utah residents (CEPH) with N. and W. European ancestry", "CHB": "Han Chinese in Beijing, China",
    "CHS": "Han Chinese South", "CLM": "Colombian in Medellin, Colombia",
    "ESN": "Esan in Nigeria", "FIN": "Finnish in Finland", "GBR": "British from England and Scotland",
    "GIH": "Gujarati Indians in Houston, TX", "GWD": "Gambian in Western Division, The Gambia",
    "IBS": "Iberian Populations in Spain", "ITU": "Indian Telugu in the UK",
    "JPT": "Japanese in Tokyo, Japan", "KHV": "Kinh in Ho Chi Minh City, Vietnam",
    "LWK": "Luhya in Webuye, Kenya", "MSL": "Mende in Sierra Leone",
    "MXL": "Mexican Ancestry in Los Angeles, CA", "PEL": "Peruvian in Lima, Peru",
    "PJL": "Punjabi in Lahore, Pakistan", "PUR": "Puerto Rican in Puerto Rico",
    "STU": "Sri Lankan Tamil in the UK", "TSI": "Toscani in Italia", "YRI": "Yoruba in Ibadan, Nigeria",
}

# Posterior threshold for "UNK", as in Martin et al. 2017 (collapse_ancestry.py).
UNK_THRESHOLD = 0.9


def group_of(population: str, super_population: str) -> str:
    """Benchmark group: AFR-American or AMR for admixed populations, else the super-population."""
    return GROUPS.get(population, super_population)


def chrom_list(value: str | None = None) -> list[int]:
    """Parse a CHROMS setting such as "1 2 22" or "1-22"."""
    value = (value or os.environ.get("CHROMS") or "1-22").replace(",", " ")
    out: list[int] = []
    for token in value.split():
        if "-" in token:
            a, b = token.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(token))
    bad = [c for c in out if c not in AUTOSOMES]
    if bad:
        raise ValueError(f"only autosomes 1-22 have pretrained models, got {bad}")
    return sorted(set(out))


def threads() -> int:
    return max(1, int(os.environ.get("THREADS", "4")))


def memory_limit_gb() -> float:
    """Memory available to this process tree: cgroup limit (Docker) or MemAvailable."""
    limits = []
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            v = Path(path).read_text().strip()
            if v.isdigit() and int(v) < 1 << 60:
                limits.append(int(v) / 1e9)
        except OSError:
            pass
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                limits.append(int(line.split()[1]) * 1024 / 1e9)
    except OSError:
        pass
    return min(limits) if limits else float("inf")


def workers_for_memory(per_worker_gb: float, reserve_gb: float = 4.0) -> int:
    """THREADS, reduced if the memory limit cannot hold that many workers."""
    n = threads()
    cap = int((memory_limit_gb() - reserve_gb) // per_worker_gb)
    if cap < n:
        log(f"  memory limit {memory_limit_gb():.0f} GB: using {max(1, cap)} workers instead of {n} "
            f"(about {per_worker_gb} GB each)")
    return max(1, min(n, cap))


def single_thread_env() -> None:
    """Pin numeric libraries to one thread. Parallelism comes from our own workers."""
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = "1"


def import_gnomix() -> None:
    """Put the pinned upstream Gnomix source on sys.path (its package is `src`)."""
    gnomix_dir = TOOLS_DIR / "gnomix"
    if not (gnomix_dir / "src" / "model.py").exists():
        raise FileNotFoundError(f"Gnomix source not found in {gnomix_dir}. Run setup_tools.sh.")
    if str(gnomix_dir) not in sys.path:
        sys.path.insert(0, str(gnomix_dir))


@contextmanager
def atomic_path(path: Path):
    """Yield a temporary path next to `path`; move it into place only on success.

    Steps skip work whose output file exists, so an output must never exist in a
    partly written state. The temporary name keeps the suffix (pandas and numpy
    choose formats from it).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".part.{path.name}")
    yield tmp
    os.replace(tmp, path)


def write_json(path: Path, obj) -> None:
    with atomic_path(path) as tmp:
        tmp.write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n")


def read_json(path: Path):
    return json.loads(Path(path).read_text())


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def inputs_stamp(chroms: list[int], use_fix) -> str:
    """Fingerprint of what karyograms and the corrected panel depend on: which
    samples use Gnofix, and the content of the calls. They are redone when it changes."""
    h = hashlib.sha256(bytes(bytearray(bool(x) for x in use_fix)))
    for c in chroms:
        h.update(f"{c}:{sha256_file(WORK / 'infer' / f'chr{c}.npz')}".encode())
    return h.hexdigest()


def stamp_is_current(directory: Path, stamp: str) -> bool:
    f = Path(directory) / ".inputs_stamp"
    return f.exists() and f.read_text().strip() == stamp


def write_stamp(directory: Path, stamp: str) -> None:
    Path(directory).mkdir(parents=True, exist_ok=True)
    (Path(directory) / ".inputs_stamp").write_text(stamp + "\n")
