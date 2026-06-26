"""
================================================================================
HCC4R 数据集入口脚本
================================================================================

功能：对 HCC4R Visium 切片执行完整的预处理 + Cell2location 反卷积流程，
      并自动调用空间生态位分析（Step 2），结果保存到 results/HCC4R/。

用法：
    python code/run_hcc4r.py

    可选：仅运行预处理（Step 1）
    python code/run_hcc4r.py --step1-only

    可选：仅运行空间分析（已有 adata_vis_post.h5ad 时跳过 Step 1）
    python code/run_hcc4r.py --step2-only

输出目录：
    results/HCC4R/                   ← Cell2location 反卷积结果
    results/HCC4R/spatial_niche/     ← 空间生态位分析结果

数据目录（Space Ranger 输出格式，与 CHC20/CHC23 结构相同）：
    data/HCC4R/                      ← HCC4R Visium 数据目录
================================================================================
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

# ── 项目根目录
_REPO_ROOT = Path(__file__).resolve().parent.parent

# ── 数据集专属配置 ────────────────────────────────────────────────────────────
SAMPLE1_NAME   = "HCC4R"
SAMPLE2_NAME   = "HCC6NR"         # 配对验证切片（可不传，设为 None 可跳过）
PATH_SCRNA     = _REPO_ROOT / "data" / "scRNA_reference.h5ad"
PATH_SAMPLE1   = _REPO_ROOT / "data" / "HCC4R"
PATH_SAMPLE2   = _REPO_ROOT / "data" / "HCC6NR"   # 若不需要配对验证，设为 None
OUTPUT_DIR     = _REPO_ROOT / "results" / SAMPLE1_NAME


def parse_args():
    p = argparse.ArgumentParser(description="HCC4R 完整分析流程（Step 1 + Step 2）")
    p.add_argument("--step1-only",   action="store_true", help="仅运行 Step 1（预处理 + 反卷积）")
    p.add_argument("--step2-only",   action="store_true", help="仅运行 Step 2（空间生态位分析）")
    p.add_argument("--no-sample2",   action="store_true", help="Step 1 中跳过 HCC6NR 配对验证切片")
    return p.parse_args()


def run_step1(with_sample2: bool = True):
    """运行 Step 1：预处理 + Cell2location 反卷积"""
    pipeline_dir = Path(__file__).parent / "pipeline"
    sys.path.insert(0, str(pipeline_dir))
    sys.path.insert(0, str(pipeline_dir.parent / "utils"))

    from run_preprocessing import main as preprocessing_main
    preprocessing_main(
        path_scrna=PATH_SCRNA,
        path_sample1=PATH_SAMPLE1,
        path_sample2=PATH_SAMPLE2 if with_sample2 else None,
        sample1_name=SAMPLE1_NAME,
        sample2_name=SAMPLE2_NAME,
        output_dir=OUTPUT_DIR,
    )


def run_step2():
    """运行 Step 2：空间生态位分析"""
    niche_script = Path(__file__).parent / "pipeline" / "run_spatial_niche_analysis.py"
    adata_path   = OUTPUT_DIR / "adata_vis_post.h5ad"
    out_dir      = OUTPUT_DIR / "spatial_niche"
    sig_out      = OUTPUT_DIR / "spatial_signature_genes.txt"

    if not adata_path.exists():
        print(f"[ERROR] adata_vis_post.h5ad not found at {adata_path}")
        print("       Please run Step 1 first: python code/run_hcc4r.py --step1-only")
        sys.exit(1)

    cmd = [
        sys.executable, str(niche_script),
        "--adata",         str(adata_path),
        "--out-dir",       str(out_dir),
        "--signature-out", str(sig_out),
    ]
    print(f"[run_hcc4r] Running Step 2: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"[ERROR] Step 2 failed with exit code {result.returncode}")
        sys.exit(result.returncode)


def main():
    args = parse_args()
    print(f"{'=' * 60}")
    print(f"  数据集: {SAMPLE1_NAME}")
    print(f"  输出目录: {OUTPUT_DIR}")
    print(f"{'=' * 60}")

    if args.step2_only:
        run_step2()
    elif args.step1_only:
        run_step1(with_sample2=not args.no_sample2)
    else:
        run_step1(with_sample2=not args.no_sample2)
        run_step2()


if __name__ == "__main__":
    main()
