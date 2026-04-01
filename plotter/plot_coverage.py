#!/usr/bin/env python3
"""
Python port of the provided Wolfram coverage plotter.
- Loads TSV coverage files matching *_hg38_ASv2.haplotagged_1000000.cov.tsv
- Normalizes per sample with a median filter
- Computes robust per-bin stats across samples
- Plots each sample with chromosome shading and outlier overlay
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def median_filter_1d(arr: np.ndarray, radius: int) -> np.ndarray:
    """Simple median filter; radius r -> window size 2r+1."""
    if radius <= 0:
        return arr.astype(float, copy=True)
    k = 2 * radius + 1
    pad = radius
    padded = np.pad(arr, pad, mode="edge")
    out = np.empty_like(arr, dtype=float)
    for i in range(len(arr)):
        window = padded[i : i + k]
        out[i] = np.median(window)
    return out


def norma(vec: np.ndarray, smoothing_radius: int = 10) -> np.ndarray:
    """Median normalize and smooth; replace zeros with 1."""
    med = np.nanmedian(vec)
    normed = vec / med if med != 0 else vec
    smoothed = median_filter_1d(normed, smoothing_radius)
    smoothed[smoothed == 0] = 1.0
    return smoothed


def compute_stats(samples: list[np.ndarray], smoothing_radius: int = 10):
    """Compute p1/p2/p1e/p2e/p3 across samples (per bin)."""
    normed = [norma(s, smoothing_radius) for s in samples]
    z = np.stack(normed, axis=1)  # bins x samples
    med = np.median(z, axis=1)
    q1 = np.quantile(z, 0.25, axis=1)
    q3 = np.quantile(z, 0.75, axis=1)
    iqr = q3 - q1
    p1 = med - 1.5 * iqr
    p2 = med + 1.5 * iqr
    p1e = med - 3.0 * iqr
    p2e = med + 3.0 * iqr
    p3 = med
    return p1, p2, p1e, p2e, p3, normed


def genome_ranges(chr_series: pd.Series):
    """Return list of (start_idx, end_idx, chr_name) for consecutive blocks."""
    ranges = []
    current_chr = chr_series.iloc[0]
    start = 0
    for i, val in enumerate(chr_series):
        if val != current_chr:
            ranges.append((start, i - 1, current_chr))
            current_chr = val
            start = i
    ranges.append((start, len(chr_series) - 1, current_chr))
    return ranges


def plot_sample(
    idx: int,
    fname: Path,
    chr_labels: pd.Series,
    p0: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    p1e: np.ndarray,
    p2e: np.ndarray,
    p3: np.ndarray,
    chr_ranges: list[tuple[int, int, str]],
    y_min: float = -1.0,
    y_max: float = 2.0,
    outdir: Path | None = None,
):
    """Plot one sample with shading, ticks, outlier overlay."""
    p3_smooth = median_filter_1d(p3, 1)
    p0_norm = p0 / p3_smooth
    p0_bad = np.where((p0_norm < p1e) | (p0_norm > p2e), p0_norm, np.nan)

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(p0_norm))

    for j, (start, end, _) in enumerate(chr_ranges):
        if j % 2 == 1:  # odd blocks shaded
            ax.axvspan(start - 0.5, end + 0.5, color="0.95", zorder=0)

    ax.fill_between(x, p1, p2, step="mid", color="0.5", alpha=0.3, linewidth=0)
    ax.step(x, p1, where="mid", color="gray", linestyle="--", alpha=0.2)
    ax.step(x, p2, where="mid", color="gray", linestyle="--", alpha=0.2)
    ax.step(x, p0_norm, where="mid", color="black", linewidth=1.5, alpha=0.7)

    ax.step(x, p0_bad, where="mid", color="red", linewidth=1.5)

    for start, end, _ in chr_ranges[:-1]:
        ax.axvline(end + 0.5, color="0.7", linestyle="--", linewidth=0.8)
    tick_positions = [(start + end) / 2 for start, end, _ in chr_ranges]
    tick_labels = [c for _, _, c in chr_ranges]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=90)

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(-0.5, len(p0_norm) - 0.5)
    ax.set_title(fname.name, fontsize=12, fontweight="bold")
    ax.set_xlabel("Genomic bins")
    ax.set_ylabel("Normalized coverage")
    fig.tight_layout()

    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        out_path = outdir / f"{fname.stem}.png"
        fig.savefig(out_path, dpi=200)
        print(f"Saved {out_path}")
    else:
        plt.show()

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Coverage plotter (Python port of Mathematica code).")
    parser.add_argument(
        "--glob",
        default="/Users/marcushansen/coverage/*_hg38_ASv2.haplotagged_1000000.cov.tsv",
        help="Glob pattern for coverage TSV files.",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output directory for PNGs (optional).")
    parser.add_argument("--smoothing-radius", type=int, default=10, help="Median filter radius (Mathematica n).")
    args = parser.parse_args()

    # glob.glob allows absolute patterns; if no wildcard, accept the single path
    if "*" in args.glob or "?" in args.glob or "[" in args.glob:
        file_paths = [Path(p) for p in sorted(glob.glob(args.glob))]
    else:
        file_paths = [Path(args.glob)]
    files = [f for f in file_paths if f.is_file()]
    if not files:
        raise SystemExit("No files found for pattern.")

    # Load all files; expect chr column then coverage columns; take last column as sample signal
    tables = [pd.read_csv(f, sep="\t", header=None) for f in files]
    chr_labels = tables[0].iloc[:, 0]
    samples = [t.iloc[:, -1].to_numpy(dtype=float) for t in tables]

    p1, p2, p1e, p2e, p3, normed = compute_stats(samples, args.smoothing_radius)
    chr_ranges = genome_ranges(chr_labels)

    for i, (f, p0) in enumerate(zip(files, normed), start=1):
        plot_sample(
            i,
            f,
            chr_labels,
            p0,
            p1,
            p2,
            p1e,
            p2e,
            p3,
            chr_ranges,
            outdir=args.output,
        )


if __name__ == "__main__":
    main()
