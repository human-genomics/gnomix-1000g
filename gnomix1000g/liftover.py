"""Point liftover of 1-based positions with a UCSC chain file (vectorised).

A chain file is a list of ungapped alignment blocks between a source assembly
(UCSC "target", here hg19) and a destination assembly (UCSC "query", here hg38).
A single base lifts only if it lies inside a block. The source sides of the
blocks in hg19ToHg38.over.chain.gz do not overlap, so each base lifts to at
most one place. We check this when the file is read.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(alleles: np.ndarray) -> np.ndarray:
    return np.array([a.translate(_COMPLEMENT) for a in alleles], dtype=alleles.dtype)


class Chain:
    def __init__(self, path: str | Path):
        rows: dict[str, list[tuple]] = {}
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            block_rows: list[tuple] = []
            for line in handle:
                f = line.split()
                if not f:
                    continue
                if f[0] == "chain":
                    if f[4] != "+":
                        raise ValueError("source strand must be '+'")
                    src_pos, dst_pos = int(f[5]), int(f[10])
                    dst_name, dst_size, dst_minus = f[7], int(f[8]), f[9] == "-"
                    block_rows = rows.setdefault(f[2], [])
                    continue
                size = int(f[0])
                block_rows.append((src_pos, src_pos + size, dst_pos, dst_minus, dst_size, dst_name))
                if len(f) == 3:
                    src_pos += size + int(f[1])
                    dst_pos += size + int(f[2])

        self._blocks = {}
        for name, b in rows.items():
            b.sort(key=lambda r: r[0])
            cols = list(zip(*b))
            src_start = np.array(cols[0], dtype=np.int64)
            src_end = np.array(cols[1], dtype=np.int64)
            if np.any(src_start[1:] < src_end[:-1]):
                raise ValueError(f"overlapping source blocks on {name}")
            self._blocks[name] = (
                src_start,
                src_end,
                np.array(cols[2], dtype=np.int64),
                np.array(cols[3], dtype=bool),
                np.array(cols[4], dtype=np.int64),
                np.array(cols[5], dtype=object),
            )

    def lift(self, chrom: str, pos) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Lift 1-based positions. Returns (chrom, pos, minus); unmapped -> ("", -1, False)."""
        pos = np.asarray(pos, dtype=np.int64)
        out_chrom = np.full(pos.shape, "", dtype=object)
        out_pos = np.full(pos.shape, -1, dtype=np.int64)
        out_minus = np.zeros(pos.shape, dtype=bool)
        if chrom not in self._blocks:
            return out_chrom, out_pos, out_minus
        src_start, src_end, dst_start, dst_minus, dst_size, dst_name = self._blocks[chrom]
        z = pos - 1
        i = np.searchsorted(src_start, z, side="right") - 1
        ic = np.clip(i, 0, None)
        hit = (i >= 0) & (z < src_end[ic])
        b = ic[hit]
        d = dst_start[b] + (z[hit] - src_start[b])
        minus = dst_minus[b]
        d = np.where(minus, dst_size[b] - 1 - d, d)
        out_chrom[hit] = dst_name[b]
        out_pos[hit] = d + 1
        out_minus[hit] = minus
        return out_chrom, out_pos, out_minus
