import argparse
from pathlib import Path

import scanpy as sc


def print_structure(adata, name: str) -> None:
    """打印单细胞或空间转录组数据的结构信息

    参数:
        adata: AnnData对象，包含单细胞或空间转录组数据
        name: str，数据集的名称，用于输出标识

    返回值:
        None
    """
    # 打印数据集标题和基本信息
    print(f"\n=== {name} ===")
    # 打印AnnData对象的摘要信息
    print(adata)
    # 打印数据矩阵的形状（样本数×特征数）
    print(f"shape: {adata.shape}")
    # 打印观察值（obs）的列名（样本元数据）
    print(f"obs columns: {list(adata.obs.columns)}")
    # 打印变量（var）的列名（特征元数据）
    print(f"var columns: {list(adata.var.columns)}")
    # 打印非结构化数据（uns）的键名
    print(f"uns keys: {list(adata.uns.keys())}")
    # 打印数据层（layers）的键名
    print(f"layers: {list(adata.layers.keys())}")



def print_query(adata_sc, adata_sp, query: str, head: int) -> None:
    if query == "sc_obs":
        print(adata_sc.obs.head(head))
    elif query == "sc_var":
        print(adata_sc.var.head(head))
    elif query == "sp_obs":
        print(adata_sp.obs.head(head))
    elif query == "sp_var":
        print(adata_sp.var.head(head))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Read scRNA + spatial data and print their structures.")
    parser.add_argument("--sc-h5ad", type=Path, required=True, help="Path to scRNA .h5ad file")
    parser.add_argument("--spatial-h5ad", type=Path, default=None, help="Path to spatial .h5ad file")
    parser.add_argument("--visium-dir", type=Path, default=None, help="Path to a Visium folder")
    parser.add_argument("--query", choices=["sc_obs", "sc_var", "sp_obs", "sp_var"], default=None)
    parser.add_argument("--head", type=int, default=5)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.spatial_h5ad is None and args.visium_dir is None:
        raise ValueError("Please provide --spatial-h5ad or --visium-dir.")

    adata_sc = sc.read_h5ad(args.sc_h5ad)
    adata_sp = sc.read_h5ad(args.spatial_h5ad) if args.spatial_h5ad else sc.read_visium(args.visium_dir)

    print_structure(adata_sc, "scRNA")
    print_structure(adata_sp, "Spatial")

    if args.query:
        print(f"\n=== query: {args.query} (head={args.head}) ===")
        print_query(adata_sc, adata_sp, args.query, args.head)


