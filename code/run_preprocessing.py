# 基于单细胞转录组（scRNA-seq）数据指导的空间转录组（Visium）数据解卷积（Deconvolution）与空间共定位分析流程
from pathlib import Path
from time import perf_counter

import logging
import sys

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch

torch.set_float32_matmul_precision('medium') # 降低浮点数的计算精度，提升训练速度

from preprocessing import (
    align_shared_genes,
    assign_treg_label,
    check_cell2location_available,
    compute_spot_cell_proportion,
    export_signatures,
    extract_cell_state_df,
    load_scrna_h5ad,
    load_visium,
    sanitize_cell_state_df,
    setup_and_train_regression_model,
    subset_t_cells,
)

ROOT_DIR = Path("/mnt/hdd/private/gyw/Code/txbb/Single_cell_train")
DATA_DIR = ROOT_DIR / "data"

PATH_SCRNA_H5AD = DATA_DIR / "scRNA_reference.h5ad"
PATH_CHC20 = DATA_DIR / "CHC20_Visium"
PATH_CHC23 = DATA_DIR / "CHC23_Visium"
OUTPUT_DIR = ROOT_DIR / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_PATH = OUTPUT_DIR / "run_preprocessing.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_PATH, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
LOGGER = logging.getLogger(__name__)


def main():
    t0 = perf_counter() # 记录流程起始时间

    regression_model_cls, cell2location_cls = check_cell2location_available() # 自定义函数，导入包内模型类
    
    # 加载数据
    adata_sc = load_scrna_h5ad(PATH_SCRNA_H5AD, t0=t0) # 自定义函数，读入单细胞数据并输出基本信息
    adata_vis = load_visium(PATH_CHC20, sample_name="CHC20") # 自定义函数，读入 Visium 数据并输出基本信息

    shared_genes, adata_sc, adata_vis_sh = align_shared_genes(adata_sc, adata_vis) # 对齐两个数据集的共同基因

    adata_t = subset_t_cells(adata_sc)
    sc.pl.dotplot(adata_t, ["CD3D", "CD4", "FOXP3", "IL2RA"], groupby="res.3", show=False) # 散点图展示 T/NK 细胞中不同簇的关键基因表达，帮助确定 Treg 簇
    plt.savefig(OUTPUT_DIR / "t_cell_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close()

    adata_sc = assign_treg_label(adata_sc, treg_clusters=("9", "9.0"))
    print(adata_sc.obs["final_celltype"].value_counts()) # 输出新的细胞类型注释统计，确认 Treg 标签已正确分配

    model = setup_and_train_regression_model(adata_sc, regression_model_cls)

    model.plot_history(50)
    history_path = OUTPUT_DIR / "regression_training_history.png"
    plt.savefig(history_path, dpi=150, bbox_inches="tight")
    plt.close()
    LOGGER.info("单细胞模型训练曲线已保存: %s", history_path)

    adata_sc = export_signatures(model, adata_sc)
    cell_state_df = sanitize_cell_state_df(extract_cell_state_df(adata_sc))

    adata_vis_sh.X = adata_vis_sh.layers["counts"].copy()
    if sp.issparse(adata_vis_sh.X):
        adata_vis_sh.X.data = np.round(adata_vis_sh.X.data)
        adata_vis_sh.X.data[adata_vis_sh.X.data < 0] = 0
    else:
        adata_vis_sh.X = np.round(adata_vis_sh.X)
        adata_vis_sh.X[adata_vis_sh.X < 0] = 0

    cell2location_cls.setup_anndata(adata_vis_sh)
    sp_model = cell2location_cls(
        adata_vis_sh,
        cell_state_df=cell_state_df,
        N_cells_per_location=30,
        detection_alpha=20,
    )
    sp_model.train(max_epochs=1000, batch_size=1000)

    sp_model.plot_history(100)
    sp_history_path = OUTPUT_DIR / "spatial_mapping_training_history.png"
    plt.savefig(sp_history_path, dpi=150, bbox_inches="tight")
    plt.close()
    LOGGER.info("空间模型训练曲线已保存: %s", sp_history_path)

    adata_vis_sh = sp_model.export_posterior(adata_vis_sh)
    # 使用后验均值（means_cell_abundance_w_sf）作为默认丰度估计，比 q05 更能反映真实丰度中心
    adata_vis_sh.obsm["cell_abundance"] = adata_vis_sh.obsm["means_cell_abundance_w_sf"]

    spot_cell_proportion = compute_spot_cell_proportion(
        adata_vis_sh, abundance_key="means_cell_abundance_w_sf"
    )
    output_csv = OUTPUT_DIR / "spot_cell_proportion.csv"
    spot_cell_proportion.to_csv(output_csv, index=False)
    LOGGER.info("Step 1 完成。Cell2location 反卷积结果已保存至: %s", output_csv)

    # 共定位标签（niche_high）由 run_spatial_niche_analysis.py 统一生成，此处不再重复定义。
    # 运行下一步：python code/run_spatial_niche_analysis.py
    LOGGER.info("Step 2（空间生态位评分）已拆分为独立脚本：code/run_spatial_niche_analysis.py")
    LOGGER.info("Step 3（可选 Wilcoxon 验证）已拆分为独立脚本：code/run_de_analysis.py")
    LOGGER.info("请依次运行上述脚本。TCGA 生存分析请在 R 中运行 code/tcga_survival_analysis.R。")

    adata_sc.write_h5ad(OUTPUT_DIR / "adata_sc_post.h5ad")
    adata_vis_sh.write_h5ad(OUTPUT_DIR / "adata_vis_post.h5ad")
    (OUTPUT_DIR / "shared_genes.txt").write_text("\n".join(shared_genes), encoding="utf-8")

    LOGGER.info("Run completed in %.2fs", perf_counter() - t0)
    LOGGER.info("Outputs saved to %s", OUTPUT_DIR)

    return adata_sc, adata_vis_sh, shared_genes


if __name__ == "__main__":
    main()
