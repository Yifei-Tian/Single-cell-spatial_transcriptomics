"""
================================================================================
脚本名称: run_joint_hcc4r_chc20.py
功能概述: HCC4R + CHC20 联合分析入口脚本 —— 完整流程（Step 1 + Step 2）
================================================================================

【联合分析背景】
    经批次效应评估（PCA/UMAP 整合可视化、细胞类型比例相关性检验），确认
    HCC4R 与 CHC20 两个数据集的 batch effect 极小，适合联合分析。

【分析流程】
    Step 1  (run_preprocessing.py → main_joint)
            用同一个 scRNA 参考训练一次 RegressionModel，得到统一 cell_state_df；
            分别对 HCC4R 和 CHC20 做 Cell2location 反卷积；
            合并两切片 spot 级结果，加 obs["sample"] 列标记来源；
            保存 adata_vis_post_joint.h5ad（联合 niche 分析输入）。

    Step 2  (run_spatial_niche_analysis.py --per-sample-neighbors)
            基于合并后的联合 AnnData 进行空间免疫抑制生态位识别；
            启用 --per-sample-neighbors 确保空间邻域只在同一切片内建立，
            不允许 HCC4R 的 spot 和 CHC20 的 spot 跨切片互为物理邻居；
            下游 cell composition、cluster、niche score、DEG 分析
            在合并后的全表（含两个样本）上进行。

【用法】
    运行完整流程（Step 1 + Step 2）：
        python code/run_joint_hcc4r_chc20.py

    仅运行 Step 1（预处理 + 反卷积 + 合并）：
        python code/run_joint_hcc4r_chc20.py --step1-only

    仅运行 Step 2（已有 adata_vis_post.h5ad 时）：
        python code/run_joint_hcc4r_chc20.py --step2-only

【输出目录】
    results/joint_HCC4R_CHC20/                - Step 1 反卷积和合并结果
    results/joint_HCC4R_CHC20/spatial_niche/  - Step 2 联合空间生态位分析结果

【数据目录（Space Ranger 输出格式）】
    data/HCC4R/      - HCC4R Visium 数据目录
    data/CHC20_Visium/ - CHC20 Visium 数据目录
================================================================================
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# ── 项目根目录（code/run_joint_hcc4r_chc20.py → 根目录）
_REPO_ROOT = Path(__file__).resolve().parent.parent

# ── 联合分析配置 ──────────────────────────────────────────────────────────────
SAMPLE_A_NAME  = "HCC4R"
SAMPLE_B_NAME  = "CHC20"
PATH_SCRNA     = _REPO_ROOT / "data" / "scRNA_reference.h5ad"
PATH_SAMPLE_A  = _REPO_ROOT / "data" / "HCC4R"
PATH_SAMPLE_B  = _REPO_ROOT / "data" / "CHC20_Visium"
OUTPUT_DIR     = _REPO_ROOT / "results" / f"joint_{SAMPLE_A_NAME}_{SAMPLE_B_NAME}"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="HCC4R + CHC20 联合分析完整流程（Step 1 + Step 2）"
    )
    p.add_argument(
        "--step1-only", action="store_true",
        help="仅运行 Step 1（预处理 + Cell2location 反卷积 + 合并）"
    )
    p.add_argument(
        "--step2-only", action="store_true",
        help="仅运行 Step 2（空间生态位联合分析）；需先完成 Step 1"
    )
    return p.parse_args()


def run_step1() -> None:
    """Step 1：联合预处理 + Cell2location 反卷积 + 合并 spot 结果"""
    pipeline_dir = Path(__file__).parent / "pipeline"
    sys.path.insert(0, str(pipeline_dir))
    sys.path.insert(0, str(pipeline_dir.parent / "utils"))

    from run_preprocessing import main_joint
    main_joint(
        path_scrna=PATH_SCRNA,
        path_sample_a=PATH_SAMPLE_A,
        path_sample_b=PATH_SAMPLE_B,
        sample_a_name=SAMPLE_A_NAME,
        sample_b_name=SAMPLE_B_NAME,
        output_dir=OUTPUT_DIR,
    )


def run_step2() -> None:
    """
    Step 2：联合空间生态位分析。

    关键参数 --per-sample-neighbors：
        确保 kNN 空间邻域只在同一切片（sample）内建立，
        不允许 HCC4R 的 spot 和 CHC20 的 spot 跨切片互为物理邻居，
        保留空间组织结构的生物学意义。
    """
    niche_script = Path(__file__).parent / "pipeline" / "run_spatial_niche_analysis.py"
    adata_path   = OUTPUT_DIR / "adata_vis_post.h5ad"
    out_dir      = OUTPUT_DIR / "spatial_niche"
    sig_out      = OUTPUT_DIR / "spatial_signature_genes.txt"

    if not adata_path.exists():
        print(f"[ERROR] adata_vis_post.h5ad not found at {adata_path}")
        print(
            "       Please run Step 1 first:\n"
            "         python code/run_joint_hcc4r_chc20.py --step1-only"
        )
        sys.exit(1)

    cmd = [
        sys.executable, str(niche_script),
        "--adata",                str(adata_path),
        "--out-dir",              str(out_dir),
        "--signature-out",        str(sig_out),
        "--per-sample-neighbors",          # 按切片内建立空间邻域（核心参数）
    ]
    print(f"[run_joint] Running Step 2: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"[ERROR] Step 2 failed with exit code {result.returncode}")
        sys.exit(result.returncode)


def main() -> None:
    args = parse_args()
    sep = "=" * 70
    print(sep)
    print(f"  联合分析: {SAMPLE_A_NAME} + {SAMPLE_B_NAME}")
    print(f"  输出目录: {OUTPUT_DIR}")
    print(sep)

    if args.step2_only:
        run_step2()
    elif args.step1_only:
        run_step1()
    else:
        run_step1()
        run_step2()


if __name__ == "__main__":
    main()