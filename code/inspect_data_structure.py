"""
================================================================================
脚本名称: inspect_data_structure.py
功能概述: 数据结构检查工具 —— 诊断与调试辅助脚本
================================================================================

【整体任务说明】
    本脚本是一个轻量级的数据诊断工具，用于在分析流程的任意阶段快速检查
    AnnData 数据对象的结构和元数据，帮助开发者验证数据格式的正确性、
    定位问题列名、确认数据维度。

    支持的操作：
      - 打印 AnnData 对象的基本信息（形状、obs/var 列名、uns 键、layers）；
      - 按需查询并展示 obs/var 表的前 N 行。

【输入文件】
    --sc-h5ad      （必填）scRNA-seq 参考数据文件，如：
                     data/scRNA_reference.h5ad
    --spatial-h5ad （二选一）已保存的空间转录组 .h5ad 文件，如：
                     data/chc20_visium.h5ad
    --visium-dir   （二选一）Space Ranger 输出目录，直接读取原始 Visium 格式，如：
                     data/CHC20_Visium/

【输出文件】
    本脚本无文件输出，所有结果均打印至标准输出（stdout）。

【调用方式】
    python inspect_data_structure.py \\
        --sc-h5ad data/scRNA_reference.h5ad \\
        --spatial-h5ad data/chc20_visium.h5ad \\
        [--query sc_obs|sc_var|sp_obs|sp_var] \\
        [--head 10]

【依赖关系】
    可在分析流程任意节点独立调用，不依赖其他自定义模块。
    建议在以下场景使用：
      - pre.py 运行后，验证 .h5ad 文件结构；
      - run_preprocessing.py 遇到列名错误时，快速定位问题。
================================================================================
"""
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


