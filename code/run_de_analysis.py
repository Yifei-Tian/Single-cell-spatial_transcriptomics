from __future__ import annotations

from pathlib import Path
import argparse
import logging
import sys

import pandas as pd
import scanpy as sc


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Step 3 DE analysis on spatial data.")
    default_root = Path(__file__).resolve().parents[1]
    default_results = default_root / "results"

    parser.add_argument(
        "--adata",
        type=Path,
        default=default_results / "adata_vis_post.h5ad",
        help="Path to adata_vis_post.h5ad",
    )
    parser.add_argument(
        "--coloc",
        type=Path,
        default=default_results / "spot_with_coloc_label.csv",
        help="Path to spot_with_coloc_label.csv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=default_results / "spatial_signature_genes.txt",
        help="Output file for signature genes",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of top genes to export",
    )
    parser.add_argument(
        "--coloc-column",
        type=str,
        default="coloc",
        help="Column name in CSV to use as coloc labels",
    )
    parser.add_argument(
        "--spot-id-column",
        type=str,
        default="spot_id",
        help="Column name in CSV to use as spot IDs",
    )
    parser.add_argument(
        "--layer",
        type=str,
        default="log1p",
        help="Layer to use for DE (falls back to X if missing)",
    )
    return parser.parse_args()


def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def _select_target_group(categories: list) -> object:
    for cand in ("1", 1, "1.0", True, "True"):
        if cand in categories:
            return cand
    return categories[-1]


def main() -> int:
    _setup_logging()
    args = _parse_args()

    logging.info("Loading adata: %s", args.adata)
    adata = sc.read_h5ad(args.adata)

    logging.info("Loading coloc labels: %s", args.coloc)
    df = pd.read_csv(args.coloc)

    if args.spot_id_column in df.columns:
        df = df.set_index(args.spot_id_column)

    if args.coloc_column not in df.columns:
        raise ValueError(f"Missing coloc column: {args.coloc_column}")

    labels = df[args.coloc_column]
    labels = labels.reindex(adata.obs_names)

    if labels.isna().any():
        missing = int(labels.isna().sum())
        logging.warning("%d spots missing coloc labels; filling with 0.", missing)
        labels = labels.fillna(0)

    adata.obs[args.coloc_column] = labels.astype("category")

    if args.layer in adata.layers:
        adata.X = adata.layers[args.layer].copy()
    else:
        logging.warning("Layer '%s' missing; DE will run on current X.", args.layer)

    sc.tl.rank_genes_groups(adata, groupby=args.coloc_column, method="wilcoxon")

    categories = list(adata.obs[args.coloc_column].cat.categories)
    if len(categories) < 2:
        raise ValueError("coloc categories < 2; cannot compute DE signature.")

    target_group = _select_target_group(categories)
    df_de = sc.get.rank_genes_groups_df(adata, group=target_group)
    genes = df_de["names"].head(args.top_n).tolist()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(genes) + "\n", encoding="utf-8")
    logging.info("Signature genes saved: %s", args.out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

