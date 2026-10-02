#!/usr/bin/env bash
# =============================================================================
# download.sh — fetch all public inputs and verify their SHA-256
# =============================================================================
# Usage:
#   bash download.sh <downloads_dir>
#
# Files (~6.4 GB in total):
#   all_hg38.pgen.zst / .pvar.zst / .psam
#       1000 Genomes high-coverage GRCh38 phased panel (3,202 samples,
#       Byrska-Bishop et al. 2022), as the official PLINK 2 resource
#       (https://www.cog-genomics.org/plink/2.0/resources). The phased genotypes
#       are the same as in the EBI 20220422 phased VCFs.
#   pretrained_gnomix_models.tar.gz
#       Pretrained Gnomix models (AI-sandbox/gnomix, Google Drive release).
#   hg19ToHg38.over.chain.gz   UCSC liftOver chain (model SNPs are GRCh37)
#   cytoBand.hg38.txt.gz       UCSC hg38 cytobands (karyogram centromeres)
#   integrated_call_samples_v3.20130502.ALL.panel
#       The 2,504 unrelated 1000 Genomes Phase 3 samples
#   {ASW,CLM,MXL,PUR}_phase1_*  1000 Genomes Phase 1 admixture working group
#       local-ancestry tracts and global proportions (comparison only)
#
# A file that exists and passes its checksum is not downloaded again.
# =============================================================================
set -euo pipefail

DEST="${1:?usage: bash download.sh <downloads_dir>}"
mkdir -p "${DEST}"

# shellcheck source=lib.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

download() {
    local name="$1" url="$2" expected="$3"
    local path="${DEST}/${name}"
    if [[ -f "${path}" && "$(sha256 "${path}")" == "${expected}" ]]; then
        echo "  [skip] ${name} (checksum ok)"
        return 0
    fi
    echo "  [download] ${name}"
    fetch "${url}" "${path}" "${expected}"
}

# ---- 1000 Genomes high-coverage GRCh38, PLINK 2 resource (Dropbox) -----------
download all_hg38.psam \
    "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1" \
    727a01b637d5e1ef02516a53069dc552b1fa46f442004ad50ea3d069c7741cdd
download all_hg38.pvar.zst \
    "https://www.dropbox.com/scl/fi/id642dpdd858uy41og8qi/all_hg38_rs_noannot.pvar.zst?rlkey=sskyiyam1bsqweujjmxqv1h55&dl=1" \
    b07449dba35891c5546c05e4182b1246c6e1adad60886b8c5b12ab806ac54d6b
download all_hg38.pgen.zst \
    "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1" \
    61112b5cda8ca9738f81c4d98aff8f6fbdc9f3f177fe6d0870077bcb5b6a8fa3

# ---- Pretrained Gnomix models (Google Drive file of AI-sandbox/gnomix) -------
# The upstream download_pretrained_models.sh uses an old Google Drive cookie
# flow that no longer works. drive.usercontent.google.com serves the same file.
download pretrained_gnomix_models.tar.gz \
    "https://drive.usercontent.google.com/download?id=1Q0zg9zqTaZUvt42uxzE_0gcfnFvjwi33&export=download&confirm=t" \
    06cc751d06ff1244dce7ffe6366e27cbedd49992f79d37b346faeed4ea748ff1

# ---- UCSC ---------------------------------------------------------------------
download hg19ToHg38.over.chain.gz \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
    5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d
download cytoBand.hg38.txt.gz \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz" \
    e514e48b0a12a516fd5231b44241a20a282f290bdd6131cadcddc6276ee26082

# ---- 1000 Genomes sample panel and Phase 1 admixture working group LAI --------
download integrated_call_samples_v3.20130502.ALL.panel \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" \
    b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858

P1="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/ancestry_deconvolution"
download README_20120604_phase1_ancestry_deconvolution "${P1}/README_20120604_phase1_ancestry_deconvolution" \
    05b9e56364781edac20d3c3865801580ed7dc8e06727e8970963ac6838a20d1c
download ASW_phase1_ancestry_deconvolution.zip "${P1}/ASW_phase1_ancestry_deconvolution.zip" \
    8745e5215cafaf41aa07c5429ca4a4040a9d8214230dfb9aff3988cbfd7a544e
download CLM_phase1_ancestry_deconvolution.zip "${P1}/CLM_phase1_ancestry_deconvolution.zip" \
    041e925f885b253177067e4f7e8253f201cbf95dbb651775bd36aa1dc7fa3b38
download MXL_phase1_ancestry_deconvolution.zip "${P1}/MXL_phase1_ancestry_deconvolution.zip" \
    d668b330f38ba1777ed7f548d7952b506b183d7de464d2f4e78f928f291b7121
download PUR_phase1_ancestry_deconvolution.zip "${P1}/PUR_phase1_ancestry_deconvolution.zip" \
    071dae616251c4dbdd255290a7579b1171622f57866596321e21ea11423653a1
download ASW_phase1_globalproportions_and_unknown.txt "${P1}/ASW_phase1_globalproportions_and_unknown.txt" \
    f503ae096a498bd8e75dbd6daf5e1449669f1fbc483a18773ee74c7533a78ae5
download CLM_phase1_globalproportions_and_unknown.txt "${P1}/CLM_phase1_globalproportions_and_unknown.txt" \
    021eaa328796304d9854da0d62d5edfa18b80fb0a2e399303555ad48b0667335
download MXL_phase1_globalproportions_and_unknown.txt "${P1}/MXL_phase1_globalproportions_and_unknown.txt" \
    b467ab56f37ef9b18ff82735e9838e4c5a7bdd6d995020b2d1738d055e4904c8
download PUR_phase1_globalproportions_and_unknown.txt "${P1}/PUR_phase1_globalproportions_and_unknown.txt" \
    c087c1464117024aa09a677e2843f06cc4a7fda81fef48cf6508e27464cb1c7a
