#!/usr/bin/env bash
# =============================================================================
# setup_tools.sh — fetch the pinned Gnomix source and the Beagle jar
# =============================================================================
# Usage:
#   bash setup_tools.sh <tools_dir>
#
# Installs into <tools_dir>:
#   gnomix/       AI-sandbox/gnomix at a pinned commit (inference code for the
#                 pretrained models, including Gnofix). Unmodified upstream code.
#   beagle.jar    Beagle 5.5 (27Feb25.75f), used only by the trio phasing
#                 benchmark. Needs Java 8 or newer.
#
# Both downloads are checked against SHA-256. Existing verified files are kept.
# =============================================================================
set -euo pipefail

DEST="${1:?usage: bash setup_tools.sh <tools_dir>}"
mkdir -p "${DEST}"

GNOMIX_COMMIT="cd15f65eadfca8e2bbf209eb69a82978a200c5fd"
GNOMIX_SHA256="846ef1ba78b569c04acaadf1c6c6e6698de881b005461f9d3ea70bf29bd6b6d8"
BEAGLE_VERSION="27Feb25.75f"
BEAGLE_SHA256="7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70"

# shellcheck source=lib.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

# ---- Gnomix ------------------------------------------------------------------
if [[ -f "${DEST}/gnomix/.commit" && "$(cat "${DEST}/gnomix/.commit")" == "${GNOMIX_COMMIT}" ]]; then
    echo "  [skip] gnomix ${GNOMIX_COMMIT:0:7}"
else
    echo "  [download] gnomix ${GNOMIX_COMMIT:0:7}"
    fetch "https://github.com/AI-sandbox/gnomix/archive/${GNOMIX_COMMIT}.tar.gz" \
        "${DEST}/gnomix.tar.gz" "${GNOMIX_SHA256}"
    staging="${DEST}/gnomix.staging"
    mkdir -p "${staging}"
    tar -xzf "${DEST}/gnomix.tar.gz" -C "${staging}"
    if [[ -d "${DEST}/gnomix" ]]; then
        mv "${DEST}/gnomix" "${DEST}/gnomix.old.$(date +%s)"
    fi
    mv "${staging}/gnomix-${GNOMIX_COMMIT}" "${DEST}/gnomix"
    rmdir "${staging}"
    echo "${GNOMIX_COMMIT}" > "${DEST}/gnomix/.commit"
fi

# ---- Beagle ------------------------------------------------------------------
if [[ -f "${DEST}/beagle.jar" && "$(sha256 "${DEST}/beagle.jar")" == "${BEAGLE_SHA256}" ]]; then
    echo "  [skip] beagle ${BEAGLE_VERSION}"
else
    echo "  [download] beagle ${BEAGLE_VERSION}"
    fetch "https://faculty.washington.edu/browning/beagle/beagle.${BEAGLE_VERSION}.jar" \
        "${DEST}/beagle.jar" "${BEAGLE_SHA256}"
fi
