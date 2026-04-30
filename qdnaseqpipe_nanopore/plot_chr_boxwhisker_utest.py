#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu


CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]
PLOT_CHR_ORDER = [str(i) for i in range(1, 23)]
BONFERRONI_TESTS = 22
ALPHA = 0.05
DELTA_MEDIAN_THRESHOLD = 0.035
ABS_MEDIAN_THRESHOLD = 0.03


def pvalue_to_stars(p_value: float) -> str:
    if pd.isna(p_value):
        return ""
    bonferroni_alpha = ALPHA / BONFERRONI_TESTS
    if p_value < bonferroni_alpha / 100:
        return "***"
    if p_value < bonferroni_alpha / 10:
        return "**"
    if p_value < bonferroni_alpha:
        return "*"
    return "ns"


def read_log2_bed(path: str, sample_name: str | None = None) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        comment="t",
        header=None,
        names=["chrom", "start", "end", "region", "log2", "strand"],
        usecols=[0, 1, 2, 3, 4, 5],
    )

    df["chrom"] = df["chrom"].astype(str).str.replace("^chr", "", regex=True)
    df = df[df["chrom"].isin(CHR_ORDER)].copy()
    df["log2"] = pd.to_numeric(df["log2"], errors="coerce")
    df = df.dropna(subset=["log2"])
    df["chrom"] = pd.Categorical(df["chrom"], categories=CHR_ORDER, ordered=True)

    if sample_name is None:
        sample_name = Path(path).stem.replace("_log2", "")
    df["sample"] = sample_name
    return df


def mann_whitney_by_chr(df1: pd.DataFrame, df2: pd.DataFrame, sample1: str, sample2: str) -> pd.DataFrame:
    results = []
    for chrom in CHR_ORDER:
        vals1 = df1.loc[df1["chrom"] == chrom, "log2"]
        vals2 = df2.loc[df2["chrom"] == chrom, "log2"]

        if len(vals1) == 0 or len(vals2) == 0:
            results.append(
                {
                    "chromosome": chrom,
                    f"n_{sample1}": len(vals1),
                    f"n_{sample2}": len(vals2),
                    "u_statistic": pd.NA,
                    "p_value": pd.NA,
                }
            )
            continue

        u_stat, p_value = mannwhitneyu(vals1, vals2, alternative="two-sided")
        results.append(
            {
                "chromosome": chrom,
                f"n_{sample1}": len(vals1),
                f"n_{sample2}": len(vals2),
                "u_statistic": u_stat,
                "p_value": p_value,
            }
        )

    return pd.DataFrame(results)


def build_summary_table(stats_df: pd.DataFrame, df1: pd.DataFrame, df2: pd.DataFrame, sample1: str, sample2: str) -> pd.DataFrame:
    med1 = df1.groupby("chrom", observed=False)["log2"].median().rename(f"median_{sample1}")
    med2 = df2.groupby("chrom", observed=False)["log2"].median().rename(f"median_{sample2}")

    summary = (
        stats_df.merge(med1, left_on="chromosome", right_index=True, how="left")
        .merge(med2, left_on="chromosome", right_index=True, how="left")
        .copy()
    )
    summary["delta_median"] = summary[f"median_{sample1}"] - summary[f"median_{sample2}"]
    summary["abs_delta_median"] = summary["delta_median"].abs()
    summary["max_abs_median"] = summary[[f"median_{sample1}", f"median_{sample2}"]].abs().max(axis=1)
    summary["bonferroni_significant"] = summary["p_value"] < (ALPHA / BONFERRONI_TESTS)
    summary["abnormal_call"] = (
        summary["bonferroni_significant"]
        & (summary["abs_delta_median"] >= DELTA_MEDIAN_THRESHOLD)
        & (summary["max_abs_median"] >= ABS_MEDIAN_THRESHOLD)
        & (summary["chromosome"].isin(PLOT_CHR_ORDER))
    )
    return summary.sort_values(["abnormal_call", "abs_delta_median", "p_value"], ascending=[False, False, True])


def plot_boxwhiskers(
    combined: pd.DataFrame,
    summary_df: pd.DataFrame,
    sample1: str,
    sample2: str,
    output_png: str,
    plot_style: str,
) -> None:
    chromosomes = [c for c in PLOT_CHR_ORDER if c in set(combined["chrom"].astype(str))]
    data1 = [combined[(combined["chrom"].astype(str) == c) & (combined["sample"] == sample1)]["log2"].tolist() for c in chromosomes]
    data2 = [combined[(combined["chrom"].astype(str) == c) & (combined["sample"] == sample2)]["log2"].tolist() for c in chromosomes]

    positions = list(range(len(chromosomes)))
    offset = 0.18
    width = 0.32

    fig, ax = plt.subplots(figsize=(max(12, len(chromosomes) * 0.6), 6))

    handles = []
    if plot_style == "box":
        bp1 = ax.boxplot(
            data1,
            positions=[p - offset for p in positions],
            widths=width,
            patch_artist=True,
            showfliers=False,
        )
        bp2 = ax.boxplot(
            data2,
            positions=[p + offset for p in positions],
            widths=width,
            patch_artist=True,
            showfliers=False,
        )

        for patch in bp1["boxes"]:
            patch.set_facecolor("#4C78A8")
            patch.set_alpha(0.8)
        for patch in bp2["boxes"]:
            patch.set_facecolor("#F58518")
            patch.set_alpha(0.8)
        handles = [bp1["boxes"][0], bp2["boxes"][0]]
    elif plot_style == "violin":
        vp1 = ax.violinplot(
            data1,
            positions=[p - offset for p in positions],
            widths=width,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )
        vp2 = ax.violinplot(
            data2,
            positions=[p + offset for p in positions],
            widths=width,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        for body in vp1["bodies"]:
            body.set_facecolor("#4C78A8")
            body.set_edgecolor("#4C78A8")
            body.set_alpha(0.5)
        for body in vp2["bodies"]:
            body.set_facecolor("#F58518")
            body.set_edgecolor("#F58518")
            body.set_alpha(0.5)
        vp1["cmedians"].set_color("#1f4e79")
        vp2["cmedians"].set_color("#a35510")
        handle1 = ax.scatter([], [], s=20, color="#4C78A8", alpha=0.8)
        handle2 = ax.scatter([], [], s=20, color="#F58518", alpha=0.8)
        handles = [handle1, handle2]
    else:
        rng = np.random.default_rng(42)
        for pos, values in enumerate(data1):
            x = rng.normal(loc=positions[pos] - offset, scale=0.04, size=len(values))
            ax.scatter(x, values, s=6, alpha=0.35, color="#4C78A8", edgecolors="none")
        for pos, values in enumerate(data2):
            x = rng.normal(loc=positions[pos] + offset, scale=0.04, size=len(values))
            ax.scatter(x, values, s=6, alpha=0.35, color="#F58518", edgecolors="none")
        median1 = [np.median(v) if len(v) else np.nan for v in data1]
        median2 = [np.median(v) if len(v) else np.nan for v in data2]
        handle1 = ax.scatter([], [], s=20, color="#4C78A8", alpha=0.8)
        handle2 = ax.scatter([], [], s=20, color="#F58518", alpha=0.8)
        ax.scatter([p - offset for p in positions], median1, s=22, color="#1f4e79", marker="_", linewidths=2)
        ax.scatter([p + offset for p in positions], median2, s=22, color="#a35510", marker="_", linewidths=2)
        handles = [handle1, handle2]

    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xticks(positions)
    ax.set_xticklabels(chromosomes)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("log2")
    bonferroni_alpha = ALPHA / BONFERRONI_TESTS
    ax.set_title(
        f"Chromosome-wise log2 distributions: {sample1} vs {sample2}\n"
        f"CALL = p < {bonferroni_alpha:.4g}, |Δmedian| ≥ {DELTA_MEDIAN_THRESHOLD}, max|median| ≥ {ABS_MEDIAN_THRESHOLD}"
    )
    ax.legend(handles, [sample1, sample2], loc="best")

    y_min = combined[combined["chrom"].astype(str).isin(chromosomes)]["log2"].min()
    y_max = combined[combined["chrom"].astype(str).isin(chromosomes)]["log2"].max()
    y_range = y_max - y_min if y_max > y_min else 1.0
    ax.set_ylim(y_min - 0.05 * y_range, y_max + 0.18 * y_range)

    stats_map = summary_df.set_index("chromosome")["p_value"].to_dict()
    abnormal_map = summary_df.set_index("chromosome")["abnormal_call"].to_dict()
    for pos, chrom in enumerate(chromosomes):
        chromosome_values = data1[pos] + data2[pos]
        if not chromosome_values:
            continue
        y_pos = max(chromosome_values) + 0.05 * y_range
        if abnormal_map.get(chrom, False):
            ax.text(pos, y_pos, "CALL", ha="center", va="bottom", fontsize=10, fontweight="bold", color="darkred")
            ax.axvspan(pos - 0.5, pos + 0.5, color="#d62728", alpha=0.08)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot chromosome-wise boxplots and run Mann-Whitney U tests.")
    parser.add_argument("bed1", help="First *_log2.bed file")
    parser.add_argument("bed2", help="Second *_log2.bed file")
    parser.add_argument("--label1", help="Label for first sample")
    parser.add_argument("--label2", help="Label for second sample")
    parser.add_argument("--out-prefix", help="Output prefix")
    parser.add_argument("--plot-style", choices=["box", "dots", "violin"], default="box", help="Plot as boxplots, jittered dots, or violins")
    args = parser.parse_args()

    sample1 = args.label1 or Path(args.bed1).stem.replace("_log2", "")
    sample2 = args.label2 or Path(args.bed2).stem.replace("_log2", "")
    out_prefix = args.out_prefix or sample1

    df1 = read_log2_bed(args.bed1, sample1)
    df2 = read_log2_bed(args.bed2, sample2)
    combined = pd.concat([df1, df2], ignore_index=True)

    output_png = f"{out_prefix}_chr_boxwhisker.png"
    output_tsv = f"{out_prefix}_mannwhitneyu.tsv"
    output_summary_tsv = f"{out_prefix}_chromosome_summary.tsv"

    stats_df = mann_whitney_by_chr(df1, df2, sample1, sample2)
    summary_df = build_summary_table(stats_df, df1, df2, sample1, sample2)
    plot_boxwhiskers(combined, summary_df, sample1, sample2, output_png, args.plot_style)
    stats_df.to_csv(output_tsv, sep="\t", index=False)
    summary_df.to_csv(output_summary_tsv, sep="\t", index=False)

    print(f"Saved plot to {output_png}")
    print(f"Saved Mann-Whitney U results to {output_tsv}")
    print(f"Saved chromosome summary to {output_summary_tsv}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
