#!/usr/bin/env python3

import argparse
import html
import math
import subprocess
from pathlib import Path

import pandas as pd


def parse_info_field(info: str) -> dict[str, str]:
    parsed = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
    return parsed


def read_log2_bed(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            rows.append(parts[:6])

    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "region", "log2", "strand"])
    if df.empty:
        return df

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["log2"] = pd.to_numeric(df["log2"], errors="coerce")
    return df.dropna(subset=["start", "end", "log2"])


def read_log2_vcf(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            info = parse_info_field(parts[7])
            rows.append(
                {
                    "chrom": parts[0],
                    "start": parts[1],
                    "end": info.get("END"),
                    "svtype": info.get("SVTYPE"),
                    "svlen": info.get("SVLEN"),
                    "bins": info.get("BINS"),
                    "score": info.get("SCORE"),
                    "log2": info.get("LOG2CNT"),
                    "filter": parts[6],
                    "alt": parts[4],
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    for column in ["start", "end", "svlen", "bins", "score", "log2"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df.dropna(subset=["start", "end", "log2"])


def read_seg(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    rename_map = {
        "SAMPLE_NAME": "sample_name",
        "CHROMOSOME": "chrom",
        "START": "start",
        "STOP": "end",
        "DATAPOINTS": "datapoints",
        "LOG2_RATIO_MEAN": "log2",
    }
    df = df.rename(columns=rename_map)
    if "sample_name" in df.columns:
        df["sample_name"] = (
            df["sample_name"]
            .astype(str)
            .str.replace(r"_segments_with_CN_class\.seg$", "", regex=True)
            .str.replace(r"_segments_log2\.seg$", "", regex=True)
        )
    for column in ["start", "end", "datapoints", "log2"]:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df.dropna(subset=[column for column in ["start", "end", "log2"] if column in df.columns])


def resolve_data_table(sample: str, input_dir: Path) -> tuple[Path, str, pd.DataFrame]:
    candidates = [
        (input_dir / f"{sample}_segments_log2.seg", "seg", read_seg),
        (input_dir / f"{sample}_segments_with_CN_class.seg", "seg", read_seg),
    ]

    for path, label, reader in candidates:
        if path.exists():
            return path, label, reader(path)

    raise SystemExit(
        "Missing SEG data table. Expected one of: "
        f"{sample}_segments_log2.seg or {sample}_segments_with_CN_class.seg"
    )


def slim_summary_df(summary_df: pd.DataFrame) -> pd.DataFrame:
    summary_df = summary_df.copy()
    if "delta_median" in summary_df.columns:
        summary_df["relative_change"] = summary_df["delta_median"].apply(lambda value: math.pow(2, value) if pd.notna(value) else pd.NA)
    if "p_value" in summary_df.columns and "abnormal_call" in summary_df.columns:
        summary_df["p_value"] = summary_df.apply(
            lambda row: row["p_value"] if str(row.get("abnormal_call", "")).lower() == "true" else "ns",
            axis=1,
        )
    preferred_columns = [
        "chromosome",
        "delta_median",
        "relative_change",
        "p_value",
        "abnormal_call",
    ]
    available = [column for column in preferred_columns if column in summary_df.columns]
    return summary_df.loc[:, available] if available else summary_df


def format_value(value, column_name: str | None = None) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if column_name == "p_value":
            return f"{value:.4e}"
        return f"{value:.4f}"
    return str(value)


def read_summary_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def dataframe_to_markdown(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    columns = [str(column) for column in df.columns]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = []
    columns = [str(column) for column in df.columns]
    for row in df.itertuples(index=False, name=None):
        values = [format_value(value, column_name) for value, column_name in zip(row, columns)]
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join([header, separator, *rows])


def dataframe_to_html(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    return df.to_html(index=False, border=0, classes="dataframe")


def build_markdown(sample: str, qdnaseq_png: Path, data_df: pd.DataFrame, data_path: Path, data_label: str, chr_box_png: Path | None, summary_df: pd.DataFrame | None, summary_path: Path | None) -> str:
    parts = [f"# Sample report: {sample}", ""]

    parts.extend([
        "## QDNAseq genome plot",
        "",
        f"![QDNAseq genome plot]({qdnaseq_png.name})",
        "",
    ])

    parts.extend([
        f"## {data_label.upper()} data table",
        "",
        f"Source: `{data_path.name}`",
        "",
        f"Rows: {len(data_df)} (showing first 100)",
        "",
        dataframe_to_markdown(data_df, max_rows=100),
        "",
    ])

    if chr_box_png is not None:
        parts.extend([
            "## Chromosome box-whisker plot",
            "",
            f"![Chromosome box-whisker plot]({chr_box_png.name})",
            "",
        ])

    if summary_df is not None and summary_path is not None:
        significant_col = None
        for candidate in ["significant", "bonferroni_significant", "abnormal_call"]:
            if candidate in summary_df.columns:
                significant_col = candidate
                break

        if significant_col is not None:
            filtered = summary_df[summary_df[significant_col].astype(str).str.lower().isin(["true", "1", "yes"])]
        else:
            filtered = summary_df

        if filtered.empty:
            filtered = summary_df
        filtered = slim_summary_df(filtered)

        parts.extend([
            f"## Significant rows from `{summary_path.name}`",
            "",
            dataframe_to_markdown(filtered),
            "",
        ])

    return "\n".join(parts)


def build_html(sample: str, qdnaseq_png: Path, data_df: pd.DataFrame, data_path: Path, data_label: str, chr_box_png: Path | None, summary_df: pd.DataFrame | None, summary_path: Path | None) -> str:
    significant_html = "<p>No chromosome summary table found.</p>"
    if summary_df is not None and summary_path is not None:
        significant_col = None
        for candidate in ["significant", "bonferroni_significant", "abnormal_call"]:
            if candidate in summary_df.columns:
                significant_col = candidate
                break

        if significant_col is not None:
            filtered = summary_df[summary_df[significant_col].astype(str).str.lower().isin(["true", "1", "yes"])]
        else:
            filtered = summary_df

        if filtered.empty:
            filtered = summary_df
        filtered = slim_summary_df(filtered)
        significant_html = f"<h2>Significant rows from <code>{html.escape(summary_path.name)}</code></h2>{dataframe_to_html(filtered)}"

    chr_box_html = ""
    if chr_box_png is not None:
        chr_box_html = (
            "<h2>Chromosome box-whisker plot</h2>"
            f'<img src="{html.escape(chr_box_png.name)}" alt="Chromosome box-whisker plot">'
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Sample report: {html.escape(sample)}</title>
  <style>
    @page {{ size: A4; margin: 2cm; }}
    body {{ font-family: Arial, sans-serif; margin: 0; }}
    h1, h2 {{ margin-bottom: 0.4em; }}
    img {{ width: 100%; max-width: 100%; height: auto; margin: 8px 0 20px; display: block; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 12px; margin: 10px 0 24px; }}
    th, td {{ border: 1px solid #ccc; padding: 4px 6px; text-align: left; }}
    th {{ background: #f5f5f5; }}
    code {{ background: #f3f3f3; padding: 1px 4px; }}
  </style>
</head>
<body>
  <h1>Sample report: {html.escape(sample)}</h1>
  <h2>QDNAseq genome plot</h2>
  <img src="{html.escape(qdnaseq_png.name)}" alt="QDNAseq genome plot">
  <h2>{html.escape(data_label.upper())} data table</h2>
  <p>Source: <code>{html.escape(data_path.name)}</code></p>
  <p>Rows: {len(data_df)} (showing first 100)</p>
  {dataframe_to_html(data_df, max_rows=100)}
  {chr_box_html}
  {significant_html}
</body>
</html>
"""


def run_pandoc(markdown_path: Path, pdf_path: Path) -> tuple[bool, str]:
    pandoc = subprocess.run(
        [
            "pandoc",
            str(markdown_path),
            "-V",
            "geometry:margin=2cm",
            "-V",
            "papersize:a4",
            "-o",
            str(pdf_path),
        ],
        capture_output=True,
        text=True,
    )
    if pandoc.returncode == 0:
        return True, ""
    return False, (pandoc.stderr or pandoc.stdout).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate markdown/pdf CNV report for a sample.")
    parser.add_argument("sample", help="Sample name, e.g. S23duffsem")
    parser.add_argument("--input-dir", default=".", help="Directory containing input files")
    parser.add_argument("--output-dir", default=".", help="Directory for report outputs")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sample = args.sample
    qdnaseq_png = input_dir / f"{sample}_QDNAseq_genome_plot.png"
    if not qdnaseq_png.exists():
        raise SystemExit(f"Missing QDNAseq plot: {qdnaseq_png}")

    data_path, data_label, data_df = resolve_data_table(sample, input_dir)

    direct_chr_box = input_dir / f"{sample}_chr_boxwhisker.png"
    if direct_chr_box.exists():
        chr_box_png = direct_chr_box
    else:
        matches = sorted(input_dir.glob(f"*{sample}*_chr_boxwhisker.png"))
        chr_box_png = matches[0] if matches else None

    summary_matches = sorted(input_dir.glob(f"*{sample}*_chromosome_summary.tsv"))
    summary_path = summary_matches[0] if summary_matches else None

    summary_df = read_summary_tsv(summary_path) if summary_path is not None else None

    md_path = output_dir / f"{sample}_report.md"
    html_path = output_dir / f"{sample}_report.html"
    pdf_path = output_dir / f"{sample}_report.pdf"

    markdown = build_markdown(sample, qdnaseq_png, data_df, data_path, data_label, chr_box_png, summary_df, summary_path)
    md_path.write_text(markdown)

    html_path.write_text(build_html(sample, qdnaseq_png, data_df, data_path, data_label, chr_box_png, summary_df, summary_path))

    try:
        pdf_ok, pdf_error = run_pandoc(md_path, pdf_path)
    except FileNotFoundError:
        pdf_ok, pdf_error = False, "pandoc not found in PATH"

    print(f"Wrote markdown: {md_path}")
    print(f"Wrote html: {html_path}")
    if pdf_ok:
        print(f"Wrote pdf: {pdf_path}")
    else:
        print(f"Skipped pdf: {pdf_error}")


if __name__ == "__main__":
    main()

