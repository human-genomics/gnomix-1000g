#!/usr/bin/env bash
# =============================================================================
# verify_pfile.sh — independent check of the phase-corrected panel
# =============================================================================
# Needs plink2 and bcftools on PATH (not used by the pipeline itself).
#
# Usage:
#   bash validation/verify_pfile.sh <data_dir> <chrom> <from_bp> <to_bp>
#   e.g. bash validation/verify_pfile.sh gnomix-1000g-data 22 20000000 22000000
#
# Exports a region from output/pfile/chrN_gnofix and from the original panel
# (downloads/all_hg38.*) and checks that samples without corrections are
# identical and that corrected samples differ only in phase.
# =============================================================================
set -euo pipefail
D="${1:?data dir}"; C="${2:?chrom}"; FROM="${3:?from bp}"; TO="${4:?to bp}"
T="$(mktemp -d)"
trap 'echo "temporary files in ${T}"' EXIT

plink2 --pfile "${D}/output/pfile/chr${C}_gnofix" vzs --chr "${C}" --from-bp "${FROM}" --to-bp "${TO}" \
    --export vcf bgz --out "${T}/new" --threads 2 >/dev/null
plink2 --pgen "${D}/work/all_hg38.pgen" --pvar "${D}/downloads/all_hg38.pvar.zst" --psam "${D}/downloads/all_hg38.psam" \
    --chr "${C}" --from-bp "${FROM}" --to-bp "${TO}" --export vcf bgz --out "${T}/old" --threads 2 >/dev/null
bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' "${T}/new.vcf.gz" > "${T}/a.txt"
bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' "${T}/old.vcf.gz" > "${T}/b.txt"
bcftools query -l "${T}/old.vcf.gz" > "${T}/samples.txt"
zcat "${D}/output/pfile/swaps.tsv.gz" | awk -F'\t' -v c="${C}" 'NR > 1 && $2 == c {print $1}' | sort -u > "${T}/corrected.txt"

python3 - "${T}" <<'PY'
import sys
t = sys.argv[1]
samples = open(f"{t}/samples.txt").read().split()
corr = set(open(f"{t}/corrected.txt").read().split())
a = [l.rstrip("\n").split("\t") for l in open(f"{t}/a.txt")]
b = [l.rstrip("\n").split("\t") for l in open(f"{t}/b.txt")]
assert len(a) == len(b) and all(x[:3] == y[:3] for x, y in zip(a, b)), "records differ"
n = {"same": 0, "other_changed": 0, "phase_only": 0, "genotype_changed": 0}
for x, y in zip(a, b):
    for s, g1, g2 in zip(samples, x[3:], y[3:]):
        if g1 == g2:
            n["same"] += 1
        elif s not in corr:
            n["other_changed"] += 1
        elif sorted(g1.split("|")) == sorted(g2.split("|")):
            n["phase_only"] += 1
        else:
            n["genotype_changed"] += 1
print(f"records {len(a)}, corrected samples {len(corr)}: {n}")
sys.exit(0 if n["other_changed"] == 0 and n["genotype_changed"] == 0 else 1)
PY
echo "PASS"
