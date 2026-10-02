# lib.sh: helpers shared by main.sh, download.sh and setup_tools.sh (sourced, not run).
# shellcheck shell=bash

if command -v sha256sum >/dev/null 2>&1; then
    sha256() { sha256sum "$1" | awk '{print $1}'; }
else
    sha256() { shasum -a 256 "$1" | awk '{print $1}'; }
fi

# fetch URL PATH SHA256: download URL to PATH and verify its checksum.
# A PATH that exists with the right checksum is kept.
fetch() {
    local url="$1" path="$2" expected="$3" actual
    if [[ -f "${path}" && "$(sha256 "${path}")" == "${expected}" ]]; then
        return 0
    fi
    curl -sSfL --retry 5 --retry-delay 10 -o "${path}.part" "${url}"
    actual="$(sha256 "${path}.part")"
    if [[ "${actual}" != "${expected}" ]]; then
        echo "ERROR: checksum mismatch for ${url}" >&2
        echo "  expected ${expected}" >&2
        echo "  got      ${actual}" >&2
        echo "  The source changed or the download failed; delete ${path}.part and retry." >&2
        exit 1
    fi
    mv "${path}.part" "${path}"
}
