"""
================================================================================
HCC4R / CHC20 跨数据集批次效应验证脚本
================================================================================

功能：
    直接复用已有 Cell2location 结果中的 spot_cell_proportion CSV，
    对 HCC4R 和 CHC20 两个数据集做跨切片对比，不再重新训练模型。

默认输入：
    results/HCC4R/spot_cell_proportion_HCC4R.csv
    results/CHC20/spot_cell_proportion_CHC20.csv

默认输出：
    results/HCC4R_CHC20_batch_validation/cross_slice_comparison/
      ├── cross_slice_mean_proportion_comparison.png
      ├── cross_slice_treg_distribution.png
      └── cross_slice_celltype_boxplot.png

用法：
    python code/run_hcc4r_chc20_batch_validation.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_HCC4R_CSV = REPO_ROOT / "results" / "HCC4R" / "spot_cell_proportion_HCC4R.csv"
DEFAULT_CHC20_CSV = REPO_ROOT / "results" / "CHC20" / "spot_cell_proportion_CHC20.csv"
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT / "results" / "HCC4R_CHC20_batch_validation" / "cross_slice_comparison"
)
DEFAULT_CELL_TYPES = ["Hepatocyte", "Treg", "T/NK", "Myeloid", "Fibroblast", "B"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Skip training and validate batch effect between HCC4R and CHC20."
    )
    parser.add_argument("--hcc4r-csv", type=Path, default=DEFAULT_HCC4R_CSV)
    parser.add_argument("--chc20-csv", type=Path, default=DEFAULT_CHC20_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--cell-types",
        nargs="+",
        default=DEFAULT_CELL_TYPES,
        help="Cell types used for cross-slice comparison.",
    )
    return parser.parse_args()


def _load_proportion_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{label} CSV not found: {path}")
    df = pd.read_csv(path)
    if "spot_id" not in df.columns:
        raise ValueError(f"{label} CSV is missing required column 'spot_id': {path}")
    return df


def _compare_slices(
    prop_hcc4r: pd.DataFrame,
    prop_chc20: pd.DataFrame,
    output_dir: Path,
    cell_types: list[str],
) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)

    available = [
        cell_type
        for cell_type in cell_types
        if cell_type in prop_hcc4r.columns and cell_type in prop_chc20.columns
    ]
    if not available:
        raise ValueError("No overlapping cell types found between HCC4R and CHC20.")

    colors = {"HCC4R": "#5b9bd5", "CHC20": "#ed7d31"}

    compare_df = pd.DataFrame(
        {
            "HCC4R": prop_hcc4r[available].mean(),
            "CHC20": prop_chc20[available].mean(),
        }
    )
    fig, ax = plt.subplots(figsize=(max(6, len(available) * 1.2), 4.5))
    compare_df.plot(kind="bar", ax=ax, color=[colors["HCC4R"], colors["CHC20"]], edgecolor="none")
    ax.set_title("Mean Cell Type Proportion: HCC4R vs CHC20")
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Mean proportion")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_dir / "cross_slice_mean_proportion_comparison.png", dpi=180)
    plt.close(fig)

    if "Treg" in available:
        fig, ax = plt.subplots(figsize=(6, 4))
        prop_hcc4r["Treg"].plot.kde(ax=ax, color=colors["HCC4R"], label="HCC4R")
        prop_chc20["Treg"].plot.kde(ax=ax, color=colors["CHC20"], label="CHC20")
        ax.set_xlabel("Treg proportion")
        ax.set_title("Treg Proportion Distribution: HCC4R vs CHC20")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "cross_slice_treg_distribution.png", dpi=180)
        plt.close(fig)

    box_types = [c for c in ["Treg", "Myeloid", "Fibroblast", "Hepatocyte"] if c in available]
    if box_types:
        records = []
        for cell_type in box_types:
            records.extend(
                {"cell_type": cell_type, "proportion": value, "slice": "HCC4R"}
                for value in prop_hcc4r[cell_type]
            )
            records.extend(
                {"cell_type": cell_type, "proportion": value, "slice": "CHC20"}
                for value in prop_chc20[cell_type]
            )

        box_df = pd.DataFrame(records)
        fig, ax = plt.subplots(figsize=(max(6, len(box_types) * 2), 4.5))
        sns.boxplot(
            data=box_df,
            x="cell_type",
            y="proportion",
            hue="slice",
            palette=colors,
            ax=ax,
            fliersize=2,
        )
        ax.set_title("Cell Type Proportion Comparison: HCC4R vs CHC20")
        ax.set_xlabel("Cell type")
        ax.set_ylabel("Proportion")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "cross_slice_celltype_boxplot.png", dpi=180)
        plt.close(fig)

    return available


def main() -> None:
    args = parse_args()
    prop_hcc4r = _load_proportion_csv(args.hcc4r_csv, "HCC4R")
    prop_chc20 = _load_proportion_csv(args.chc20_csv, "CHC20")
    available = _compare_slices(
        prop_hcc4r=prop_hcc4r,
        prop_chc20=prop_chc20,
        output_dir=args.output_dir,
        cell_types=args.cell_types,
    )

    print("=" * 72)
    print("HCC4R / CHC20 batch-effect validation completed")
    print(f"HCC4R CSV : {args.hcc4r_csv}")
    print(f"CHC20 CSV : {args.chc20_csv}")
    print(f"Output dir: {args.output_dir}")
    print(f"Cell types used: {', '.join(available)}")
    print("Generated files:")
    print(f"  - {args.output_dir / 'cross_slice_mean_proportion_comparison.png'}")
    print(f"  - {args.output_dir / 'cross_slice_treg_distribution.png'}")
    print(f"  - {args.output_dir / 'cross_slice_celltype_boxplot.png'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
