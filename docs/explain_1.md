# `run_hcc4r.py` 输出文件完整说明

本文档列举运行 `python code/run_hcc4r.py` 后生成的所有结果文件（默认完整流程：Step 1 + Step 2），并附带详细解释。若额外传入 `--paper-figures` 或 `--figures-only`，还会生成 Step 3 论文图表。

---

## 一、目录结构总览

```
results/
├── HCC4R/                            ← Step 1 主输出目录
│   ├── adata_vis_post.h5ad
│   ├── adata_vis_post_HCC4R.h5ad
│   ├── adata_sc_post.h5ad
│   ├── adata_sc_post_HCC4R.h5ad
│   ├── spot_cell_proportion_HCC4R.csv
│   ├── shared_genes_HCC4R.txt
│   ├── regression_training_history_HCC4R.png
│   ├── spatial_mapping_training_history_HCC4R.png
│   ├── t_cell_dotplot_horizontal.png
│   ├── scrna_tsne_celltype.png
│   ├── scrna_celltype_marker_heatmap.png
│   ├── run_preprocessing.log
│   │
│   ├── cross_slice_comparison/       ← 有 HCC1R 验证切片时生成
│   │   ├── cross_slice_mean_proportion_comparison.png
│   │   ├── cross_slice_treg_distribution.png
│   │   └── cross_slice_celltype_boxplot.png
│   │
│   └── spatial_niche/                ← Step 2 主输出目录
│       ├── spatial_niche_scores.csv
│       ├── spatial_niche_parameters.csv
│       ├── sensitivity_analysis.csv
│       ├── neighborhood_cluster_stats.csv
│       ├── immunosuppressive_niche_signature_genes_ranked.csv
│       ├── immunosuppressive_niche_signature_genes.txt
│       ├── prior_gene_set_auc.csv
│       ├── gini_score_genes.csv
│       ├── param_scan_deg_stability.csv
│       └── plots/
│           ├── spatial_hepatocyte.png
│           ├── spatial_treg.png
│           ├── spatial_myeloid.png
│           ├── spatial_fibroblast.png
│           ├── spatial_immunosuppressive_niche_score.png
│           ├── spatial_neighborhood_clusters.png
│           ├── spatial_niche_semantic_labels.png
│           ├── spatial_region_labels.png
│           ├── spatial_niche_high_score_spots.png
│           ├── spatial_niche_cluster_vs_score_comparison.png
│           ├── distance_to_hep_high_vs_niche_score.png
│           ├── region_score_boxplots.png
│           ├── hepatocyte_vs_treg_niche_score.png
│           ├── celltype_niche_correlation.png
│           ├── lr_communication_heatmap.png
│           ├── sensitivity_niche_stability.png
│           ├── niche_signature_volcano.png
│           ├── niche_fraction_scatter.png
│           └── param_scan_deg_stability_heatmap.png
│
├── HCC4R/spatial_signature_genes.txt ← Step 2 TCGA 接口文件
│
└── paper_figures/                    ← Step 3（需 --paper-figures）
    ├── fig1A_workflow_diagram.png
    ├── fig1B_scrna_tsne_celltype.png
    ├── fig1C_scrna_marker_heatmap.png
    ├── fig1D_treg_dotplot.png
    ├── fig2A_HE_image.png
    ├── fig2B_spatial_Malignant.png
    ├── fig2B_spatial_HSC.png
    ├── fig2B_spatial_Treg.png
    ├── fig2B_spatial_Myeloid.png
    ├── fig2C_coloc_heatmap.png
    ├── fig3A_spatial_domain_map.png
    ├── fig3B_domain_composition_barplot.png
    ├── fig3C_violin_treg_like_score.png
    ├── fig3C_violin_immune_stroma_score.png
    ├── fig3C_violin_immunosuppressive_niche_score.png
    ├── fig4A_volcano.png
    ├── fig4B_niche_marker_heatmap.png
    ├── fig4C_pathway_bubble.png
    ├── fig5A_label_transfer.png
    ├── fig5B_validation_violin_*.png
    ├── fig5C_sensitivity_stability.png
    └── fig5C_param_scan_heatmap.png
```

---

## 二、Step 1 输出文件详解（`results/HCC4R/`）

> **来源脚本**：`code/pipeline/run_preprocessing.py`，由 `run_hcc4r.py` 自动调用。
> **Step 1 完成反卷积前的核心任务**：加载 scRNA 参考 → Treg 鉴定注释 → RegressionModel 训练 → Cell2location 空间建模。

---

### 数据文件

#### `adata_vis_post.h5ad`

| 属性 | 说明 |
|------|------|
| **格式** | AnnData（HDF5） |
| **内容** | HCC4R Visium 切片经 Cell2location 反卷积后的完整空间转录组对象 |
| **主要字段** | `obs`：spot 注释；`obsm["spatial"]`：空间坐标；`obsm["means_cell_abundance_w_sf"]`：Cell2location 后验细胞丰度矩阵（每个 spot 各细胞类型的估计细胞数量）；`obsm["cell_abundance"]`：同上的副本；`X`：原始 count 矩阵 |
| **用途** | Step 2 空间生态位分析的核心输入文件，是最重要的中间结果 |

#### `adata_vis_post_HCC4R.h5ad`

内容与 `adata_vis_post.h5ad` 完全相同，区别在于文件名带有样本名称后缀，是带名称标识的副本，方便多数据集并行时按名称区分。

#### `adata_sc_post.h5ad`

| 属性 | 说明 |
|------|------|
| **格式** | AnnData（HDF5） |
| **内容** | scRNA-seq 参考数据，包含 RegressionModel 训练后导出的细胞类型签名（后验参数） |
| **主要字段** | `obs["final_celltype"]`：Treg 重注释后的最终细胞类型标签；`obsm["X_tsne"]` 或 `obsm["X_umap"]`：降维坐标 |
| **用途** | 供 Step 3 论文图表脚本读取，生成 scRNA tSNE 图和 Marker 热图 |

#### `adata_sc_post_HCC4R.h5ad`

内容与 `adata_sc_post.h5ad` 完全相同，带样本名称后缀的副本。

#### `spot_cell_proportion_HCC4R.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV 表格 |
| **内容** | HCC4R 每个 spot 的各细胞类型归一化比例（行=spot，列=细胞类型，值∈[0,1]，行和为 1） |
| **列** | `spot_id`、`Hepatocyte`、`Malignant`、`Treg`、`T/NK`、`Myeloid`、`Fibroblast`、`HSC`、`B cell`、`Plasma cell` 等 |
| **用途** | 多切片对比图（cross_slice_comparison）的输入数据；供研究者直接查看每个 spot 的细胞组成 |

#### `shared_genes_HCC4R.txt`

| 属性 | 说明 |
|------|------|
| **格式** | 纯文本，每行一个基因名 |
| **内容** | scRNA-seq 参考数据与 HCC4R Visium 数据的共享基因交集 |
| **用途** | 记录实际参与 Cell2location 建模的基因集合；可用于检查共享基因数量是否合理 |

---

### 图表文件

#### `regression_training_history_HCC4R.png`

| 属性 | 说明 |
|------|------|
| **内容** | scRNA RegressionModel（参考签名学习阶段）的训练损失曲线 |
| **解读** | 横轴为 epoch，纵轴为训练损失；曲线平滑下降并趋于收敛说明训练正常 |
| **用途** | 确认 RegressionModel 训练质量，训练发散或不收敛时需调整超参数 |

#### `spatial_mapping_training_history_HCC4R.png`

| 属性 | 说明 |
|------|------|
| **内容** | Cell2location 空间建模阶段（反卷积）的训练损失曲线 |
| **解读** | 最后 100 个 epoch 的损失曲线；前 100 个 epoch 省略以突出收敛阶段 |
| **用途** | 验证 Cell2location 模型收敛，是反卷积结果可信度的直接证据 |

#### `t_cell_dotplot_horizontal.png`

| 属性 | 说明 |
|------|------|
| **内容** | T/NK 细胞亚群中 Treg 标志基因（CD3D、CD4、FOXP3、IL2RA）的横版气泡图 |
| **横轴** | 细胞聚类簇（res.3 分辨率下的 Leiden 聚类） |
| **纵轴** | Treg 标志基因 |
| **气泡大小** | 各聚类中表达该基因的细胞比例（检出率） |
| **气泡颜色** | 平均表达量（Blues 色板，颜色越深表达越高） |
| **用途** | 辅助确认哪些 T/NK 聚类为 Treg 亚群（高 FOXP3 + 高 IL2RA 的簇即被标注为 Treg）；对应论文 Figure 1D |

#### `scrna_tsne_celltype.png`

| 属性 | 说明 |
|------|------|
| **内容** | scRNA-seq 数据的 tSNE（或 UMAP）降维图，每个细胞按类型着色 |
| **风格** | 高饱和度离散色板，每种细胞类型在图中心标注名称，右侧图例；带坐标轴箭头（仿论文风格） |
| **颜色** | 每种细胞类型固定颜色（T cell=绿、Myeloid=紫、Malignant=蓝等） |
| **用途** | 展示 scRNA 参考数据的细胞类型分布格局；对应论文 Figure 1B |

#### `scrna_celltype_marker_heatmap.png`

| 属性 | 说明 |
|------|------|
| **内容** | 各细胞类型 Marker 基因的平均表达热图 |
| **纵轴** | Marker 基因（按细胞类型分组排列：T cell、Myeloid、NK、B cell、Malignant、Endothelial、Epithelial、Plasma cell、HSC） |
| **横轴** | 细胞类型 |
| **颜色** | Z-score 标准化后的平均表达量，使用 Yellow-Black-Purple 三色渐变（亮黄=高，黑色=中，深紫=低） |
| **用途** | 验证各细胞类型的 Marker 基因特异性表达；对应论文 Figure 1C |

#### `run_preprocessing.log`

| 属性 | 说明 |
|------|------|
| **内容** | Step 1 全流程运行日志，包含各步骤的时间戳、共享基因数、细胞类型计数、警告信息等 |
| **用途** | 排查运行错误；记录实际使用的参数（共享基因数、训练 epoch 等） |

---

### 多切片验证图（`cross_slice_comparison/`，有 HCC1R 数据时生成）

#### `cross_slice_mean_proportion_comparison.png`

| 属性 | 说明 |
|------|------|
| **内容** | HCC4R 与 HCC1R 各细胞类型均值比例的并排条形图 |
| **横轴** | 细胞类型（Hepatocyte、Treg、T/NK、Myeloid、Fibroblast、B cell） |
| **纵轴** | 均值比例 |
| **颜色** | 蓝色=HCC4R，橙色=HCC1R |
| **用途** | 评估两切片的细胞组成是否相似，是 batch effect 评估的核心依据 |

#### `cross_slice_treg_distribution.png`

| 属性 | 说明 |
|------|------|
| **内容** | HCC4R 与 HCC1R 的 Treg 比例分布的核密度估计（KDE）对比图 |
| **横轴** | Treg proportion |
| **用途** | 比较两切片 Treg 浸润水平的分布形状；分布相似则说明两切片免疫微环境一致 |

#### `cross_slice_celltype_boxplot.png`

| 属性 | 说明 |
|------|------|
| **内容** | Treg、Myeloid、Fibroblast、Hepatocyte 四种细胞类型比例的箱线图对比（HCC4R vs HCC1R） |
| **用途** | 量化比较两切片核心细胞类型的比例差异，中位数和四分位距接近时可认为 batch effect 极小 |

---

## 三、Step 2 输出文件详解（`results/HCC4R/spatial_niche/`）

> **来源脚本**：`code/pipeline/run_spatial_niche_analysis.py`，由 `run_hcc4r.py` 通过 subprocess 调用。
> **Step 2 核心任务**：空间邻域 kNN 图 → Leiden 聚类 → 多层次免疫抑制评分 → 签名基因提取 → 参数扫描。

---

### CSV 数据文件

#### `spatial_niche_scores.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，行=每个 Visium spot |
| **列说明** | |
| `spot_id` | spot 唯一标识符 |
| `spatial_x`, `spatial_y` | spot 的空间坐标 |
| `Hepatocyte`, `Treg`, `T/NK`, `Myeloid`, `Fibroblast`, ... | 各细胞类型的归一化比例 |
| `immunosuppressive_gene_score` | 免疫抑制基因模块评分：FOXP3/IL2RA/CTLA4/TIGIT/TGFB1 等先验基因在该 spot 的平均 CP10K+log1p 表达量 |
| `Treg_like_score` | Treg-like 评分 = Treg 比例 Z-score + immunosuppressive_gene_score Z-score |
| `immune_stroma_score` | 免疫-基质评分 = Treg + T/NK + Myeloid + Fibroblast 比例 Z-score 加总 |
| `immunosuppressive_niche_score` | **综合免疫抑制生态位评分（核心指标）** = Treg Z-score + Myeloid Z-score + Fibroblast Z-score + gene_score Z-score（五维 Z-score 加总） |
| `niche_high` | 布尔标签：niche score ≥ 第 85 百分位数阈值为 True（即最高免疫抑制的 15% spot） |
| `spatial_region` | 辅助空间区域标注（`tumor_core`/`stroma_immune`/`tumor_edge`/`other`），基于 Hepatocyte 高/低和 niche_high 规则推断 |
| `neighborhood_cluster` | Leiden 无监督聚类簇编号（基于 kNN 空间邻域组成向量聚类） |
| `niche_semantic_label` | 语义注释标签（如 `immunosuppressive_niche`/`hepatocyte_region`/`immune_stromal` 等），由 Leiden 聚类簇的细胞组成特征自动推断 |
| `distance_to_hep_high` | 该 spot 到最近高肝细胞密度 spot 的欧氏距离（辅助刻画肿瘤边界） |
| **用途** | 最核心的结果表，后续所有可视化和统计分析均基于此表；也是 `run_paper_figures.py` 的主要输入 |

#### `spatial_niche_parameters.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，单行参数记录 |
| **内容** | 本次分析使用的所有关键参数：`n_neighbors=15`、`niche_high_quantile=0.85`、`leiden_resolution=0.5`、`radius_multiplier=1.25`、`n_spots` 等 |
| **用途** | 记录分析参数供复现；确保后续比较时参数一致 |

#### `sensitivity_analysis.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，每行为一种参数扰动配置 |
| **主要列** | `param_type`（k 或 radius）、`param_value`（扰动值）、`spearman_rho`（与参考配置的 Spearman 秩相关系数）、`weighted_jaccard`（与参考配置的加权 Jaccard 相似度） |
| **解读** | `spearman_rho` 越接近 1，说明不同 k 值下 niche 评分排序越一致；`weighted_jaccard` 越高说明高分区域重叠越多，结果越稳健 |
| **用途** | 量化评估分析结果对 kNN 邻居数和半径参数的敏感性，是方法学稳健性的核心证据 |

#### `neighborhood_cluster_stats.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，行=Leiden 聚类簇 |
| **内容** | 每个邻域聚类簇的统计信息：spot 数量（n_spots）、各细胞类型的均值比例、平均 niche score 等 |
| **用途** | 理解各 Leiden 聚类簇的生物学特征；辅助人工核查聚类结果是否合理 |

#### `immunosuppressive_niche_signature_genes_ranked.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，行=基因（按综合评分降序排列） |
| **主要列** | |
| `gene` | 基因名 |
| `mean_high` | niche_high 组的平均 CP10K+log1p 表达量 |
| `mean_low` | niche_low 组的平均 CP10K+log1p 表达量 |
| `log2_fc` | 两组之间的 log2 折叠变化（niche_high / niche_low） |
| `frac_high` | niche_high 组中检出该基因（表达>0）的 spot 比例 |
| `frac_low` | niche_low 组中检出该基因（表达>0）的 spot 比例 |
| `delta_frac` | 检出率差值 = frac_high - frac_low（刻画基因在 niche_high 中的富集程度） |
| `composite_score` | 综合排序分 = 0.5×norm(log2FC) + 0.5×norm(delta_frac)（同时考虑表达量差异和检出率差异） |
| `pvalue` | Wilcoxon 秩和检验 p 值 |
| `fdr` | Benjamini-Hochberg 多重校正 FDR |
| **筛选标准** | FDR < 0.05 且（log2FC > 0.5 OR delta_frac > 0.10） |
| **用途** | 最终免疫抑制 niche 特征基因的完整排名表，是方法学分析的核心输出；供 TCGA ssGSEA 分析、下游通路富集分析使用 |

#### `immunosuppressive_niche_signature_genes.txt`

| 属性 | 说明 |
|------|------|
| **格式** | 纯文本，每行一个基因名 |
| **内容** | 从 `ranked.csv` 中提取的 Top-N（默认 80 个）特征基因列表 |
| **用途** | 直接作为 TCGA ssGSEA 预后分析的基因集输入（`tcga_survival_analysis.R` 读取此文件） |

#### `spatial_signature_genes.txt`

| 属性 | 说明 |
|------|------|
| **路径** | `results/HCC4R/spatial_signature_genes.txt`（在 spatial_niche 目录外一级） |
| **格式** | 纯文本，每行一个基因名 |
| **内容** | 同 `immunosuppressive_niche_signature_genes.txt`，是通过 `--signature-out` 参数单独指定路径的副本 |
| **用途** | 标准化接口文件，专供 TCGA 生存分析脚本读取 |

#### `prior_gene_set_auc.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，行=先验功能基因集中的基因 |
| **主要列** | `gene`、`gene_set`（Treg_markers/TAM_features/CAF_activation）、`auc`（Mann-Whitney AUC，衡量该基因在 niche_high 中的富集程度）、`frac_high`、`frac_low`、`delta_frac` |
| **解读** | AUC > 0.6 且 delta_frac > 0.05 的先验基因被强制纳入签名基因集（即使 Wilcoxon FDR 未达阈值） |
| **用途** | 记录 Layer 2 先验基因集检验结果，解释为什么 FOXP3 等低表达基因也能被纳入签名 |

#### `gini_score_genes.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，行=基因（Gini Index 筛选通过的基因） |
| **主要列** | `gene`、`gini`（Gini 不均匀性指数，衡量基因是否局灶性高表达）、`log2_fc`、`mean_high` |
| **解读** | Gini > 0.3 且 log2FC > 0 的基因被认为在少数 spot 中高度富集（局灶性表达）；FOXP3 等稀有 Treg 基因典型 Gini ≈ 0.3~0.5 |
| **用途** | Layer 3 特异性基因筛选结果；即使基因整体均值很低，只要在部分 niche_high spot 中局灶性高表达，也能被发现 |
| **注意** | 若无基因通过筛选，文件内容为空但文件仍会生成（保证路径可预期） |

#### `param_scan_deg_stability.csv`

| 属性 | 说明 |
|------|------|
| **格式** | CSV，每行为一种参数组合（k × quantile，共 5×4=20 行） |
| **主要列** | `k`（kNN 邻居数）、`quantile`（niche_high 分位数阈值）、`n_sig_deg`（FDR<0.05 且 log2FC>0.5 的 DEG 数量）、`mean_log2fc_topN`（Top-N 基因的平均 log2FC）、`top_genes_str`（Top-20 基因名，分号分隔） |
| **用途** | 参数扫描热图的数据来源；辅助选择最优 k 和 quantile 参数组合（选颜色最深的高原区域） |

---

### 可视化图表（`plots/`）

#### 空间分布图（单细胞类型）

| 文件名 | 内容 | 解读 |
|--------|------|------|
| `spatial_hepatocyte.png` | Hepatocyte 比例的空间散点图 | 颜色越深（亮黄）代表该 spot 肝细胞比例越高，反映肿瘤实质区域分布 |
| `spatial_treg.png` | Treg 比例的空间散点图 | 高 Treg 区域（亮黄）提示局部免疫抑制微环境 |
| `spatial_myeloid.png` | Myeloid 比例的空间散点图 | Myeloid 细胞富集区域与 Treg 高度重叠时提示 TAM 参与免疫抑制 |
| `spatial_fibroblast.png` | Fibroblast 比例的空间散点图 | 高成纤维细胞区域反映基质富集，是 CAF 参与免疫抑制的空间线索 |

> 以上图均使用 Yellow-Black-Purple（YBP）三色渐变：亮黄=高值，黑色=中间值，深紫=低值

#### `spatial_immunosuppressive_niche_score.png`

| 属性 | 说明 |
|------|------|
| **内容** | `immunosuppressive_niche_score` 的空间分布图（连续值着色） |
| **颜色** | YBP 渐变，亮黄区域为最高免疫抑制 spot |
| **意义** | 直观展示免疫抑制 niche 在肿瘤组织中的空间位置和密度；是整个分析的核心可视化 |

#### `spatial_neighborhood_clusters.png`

| 属性 | 说明 |
|------|------|
| **内容** | Leiden 无监督聚类结果的空间分布图（每个 spot 按聚类簇着色） |
| **意义** | 展示基于 kNN 空间邻域组成向量聚类识别出的 cellular neighborhoods，是 niche 发现的无监督主体结果 |

#### `spatial_niche_semantic_labels.png`

| 属性 | 说明 |
|------|------|
| **内容** | 对 Leiden 聚类簇进行语义注释后的空间分布图 |
| **标签类型** | `immunosuppressive_niche`（免疫抑制生态位）、`hepatocyte_region`（肝细胞富集区）、`immune_stromal`（免疫基质混合区）等 |
| **意义** | 将无监督聚类结果转化为有生物学含义的空间域标注 |

#### `spatial_region_labels.png`

| 属性 | 说明 |
|------|------|
| **内容** | 辅助空间区域标注图（规则化分区） |
| **分区** | `tumor_core`（肿瘤核心：高 Hepatocyte）、`stroma_immune`（基质-免疫区：低 Hep + 高 niche score）、`tumor_edge`（肿瘤边缘：中间 Hep）、`other` |
| **意义** | 辅助可视化，仅用于区域评分对比和箱线图分析，不参与 niche 主体发现 |

#### `spatial_niche_high_score_spots.png`

| 属性 | 说明 |
|------|------|
| **内容** | `niche_high=True`（niche score ≥ 85th percentile）的 spot 空间分布二值图 |
| **颜色** | 高免疫抑制 spot 与背景对比显示 |
| **意义** | 展示高免疫抑制 niche 的空间位置，是签名基因提取分组的可视化验证 |

#### `spatial_niche_cluster_vs_score_comparison.png`

| 属性 | 说明 |
|------|------|
| **内容** | 并排对比图：左图=Leiden 聚类法识别的 niche，右图=评分阈值法（niche_high）识别的 niche，并标注 Jaccard 相似度 |
| **意义** | 双重验证：两种独立方法（无监督聚类 vs 规则化阈值）高度重叠（高 Jaccard）说明 niche 发现结果稳健 |

#### `distance_to_hep_high_vs_niche_score.png`

| 属性 | 说明 |
|------|------|
| **内容** | 折线图：X 轴=到最近高肝细胞 spot 的距离（分组），Y 轴=平均 immunosuppressive_niche_score |
| **预期模式** | 距离肿瘤实质越近，免疫抑制评分越高（T-I 边界富集免疫抑制） |
| **意义** | 描述 niche 与肿瘤实质的空间关系，揭示免疫抑制 niche 的空间分布规律 |

#### `region_score_boxplots.png`

| 属性 | 说明 |
|------|------|
| **内容** | 各空间区域（tumor_core/stroma_immune/tumor_edge/other）的三种评分箱线图 |
| **意义** | 定量比较不同空间区域的免疫抑制程度差异；`stroma_immune` 区预期评分最高 |

#### `hepatocyte_vs_treg_niche_score.png`

| 属性 | 说明 |
|------|------|
| **内容** | 散点图：X 轴=Hepatocyte 比例，Y 轴=Treg 比例，颜色=niche_score |
| **意义** | 展示肿瘤细胞、Treg 和 niche score 的三方关系；高 Treg + 中等 Hep 的区域通常是 niche 核心 |

#### `celltype_niche_correlation.png`

| 属性 | 说明 |
|------|------|
| **内容** | 细胞类型比例与 niche score 的 Spearman 相关性热图 |
| **意义** | 量化哪些细胞类型与免疫抑制 niche 评分最相关；Treg、Myeloid 预期有最强正相关 |

#### `lr_communication_heatmap.png`

| 属性 | 说明 |
|------|------|
| **内容** | 预设配体-受体通讯对（CCL22-CCR4、TGFB1-TGFBR1、CXCL12-CXCR4 等 12 对）在 niche_high vs niche_low 的通讯强度热图 |
| **颜色** | 每格颜色=niche_high 中的 product score 均值，行=配体，列=受体 |
| **意义** | 探索免疫抑制 niche 中活跃的细胞间通讯信号轴 |

#### `sensitivity_niche_stability.png`

| 属性 | 说明 |
|------|------|
| **内容** | 2×2 的敏感性分析稳健性图：左列=Spearman ρ，右列=加权 Jaccard；上行=k 参数扰动（k=10/15/20），下行=半径倍增系数扰动（1.0/1.25/1.5） |
| **解读** | Spearman ρ > 0.8 说明不同参数下 niche 评分排序高度一致；加权 Jaccard > 0.6 说明高分区域重叠充分 |
| **意义** | 方法学稳健性的关键证据，说明 niche 识别结果不是参数敏感的偶然产物 |

#### `niche_signature_volcano.png`

| 属性 | 说明 |
|------|------|
| **内容** | 全基因火山图：X 轴=log2(FC)（niche_high vs niche_low），Y 轴=-log10(FDR) |
| **颜色** | 红色=显著上调（FC>0.5 且 FDR<0.05），蓝色=显著下调，灰色=不显著；FOXP3/TGFB1/FAP 等先验基因标注名称 |
| **意义** | 总览签名基因的统计显著性分布；火山图右上角基因为 niche 特异性表达最强的基因 |

#### `niche_fraction_scatter.png`

| 属性 | 说明 |
|------|------|
| **内容** | 检出率差值散点图：X 轴=delta_frac（niche_high 检出率 - niche_low 检出率），Y 轴=log2FC |
| **意义** | 补充火山图的不足——FOXP3 等因细胞稀释效应导致均值 log2FC 接近 0，但 delta_frac 较高，此图能将其突出显示；右上角基因是 niche 富集最可靠的候选基因 |

#### `param_scan_deg_stability_heatmap.png`

| 属性 | 说明 |
|------|------|
| **内容** | 参数扫描热图：横轴=niche_high_quantile（0.80/0.85/0.90/0.95），纵轴=k（8/10/15/20/25），颜色=FDR<0.05 且 log2FC>0.5 的 DEG 数量 |
| **颜色** | 黄→橙→红，颜色越深 DEG 越多 |
| **选参标准** | 颜色最深且处于"高原"（相邻格子数值相近）的参数组合为推荐；k=15, quantile=0.90 通常为最优区域 |
| **注意** | 0.95 列虽然最深但不是高原，是悬崖边缘（样本量过小导致假阳性），不宜选取 |

---

## 四、Step 3 论文图表（`results/paper_figures/`，需 `--paper-figures`）

> **来源脚本**：`code/pipeline/run_paper_figures.py`，调用 `paper_plot_functions.py` 中的绘图函数。

---

### Figure 1：研究概览

| 文件名 | 对应图 | 内容 |
|--------|--------|------|
| `fig1A_workflow_diagram.png` | Figure 1A | 技术路线图（纯 matplotlib 绘制，无需数据文件），展示从 Visium 空间转录组到 niche 识别再到 TCGA 预后的完整流程 |
| `fig1B_scrna_tsne_celltype.png` | Figure 1B | 复制自 `results/HCC4R/scrna_tsne_celltype.png`，scRNA-seq tSNE 细胞类型分布图 |
| `fig1C_scrna_marker_heatmap.png` | Figure 1C | 复制自 `results/HCC4R/scrna_celltype_marker_heatmap.png`，细胞类型 Marker 基因热图 |
| `fig1D_treg_dotplot.png` | Figure 1D | 复制自 `results/HCC4R/t_cell_dotplot_horizontal.png`，Treg 标志基因气泡图 |

---

### Figure 2：空间细胞组成

| 文件名 | 对应图 | 内容 |
|--------|--------|------|
| `fig2A_HE_image.png` | Figure 2A | HCC4R 切片 H&E 染色图（需提供 `--hcc4r-visium-dir` 原始目录；否则用占位图替代）；旁附 niche_high 二值叠加 |
| `fig2B_spatial_Malignant.png` | Figure 2B | Malignant 细胞空间丰度散点图（论文版，带色标） |
| `fig2B_spatial_HSC.png` | Figure 2B | HSC（肝星状细胞）空间丰度图 |
| `fig2B_spatial_Treg.png` | Figure 2B | Treg 空间丰度图 |
| `fig2B_spatial_Myeloid.png` | Figure 2B | Myeloid 细胞空间丰度图 |
| `fig2C_coloc_heatmap.png` | Figure 2C | 细胞类型共定位 Pearson 相关性热图；颜色越红说明两种细胞类型空间共定位程度越高 |

---

### Figure 3：空间域与评分

| 文件名 | 对应图 | 内容 |
|--------|--------|------|
| `fig3A_spatial_domain_map.png` | Figure 3A | 空间域地图（tumor_core/stroma_immune/tumor_edge/other 四类域着色），高分辨率论文版 |
| `fig3B_domain_composition_barplot.png` | Figure 3B | 各空间域细胞类型组成堆叠条形图（或百分比条形图），展示各域的细胞组成差异 |
| `fig3C_violin_treg_like_score.png` | Figure 3C | Treg-like Score 在各空间域的小提琴图（内嵌箱线），含显著性标注（tumor_core vs tumor_edge） |
| `fig3C_violin_immune_stroma_score.png` | Figure 3C | Immune-Stroma Score 小提琴图 |
| `fig3C_violin_immunosuppressive_niche_score.png` | Figure 3C | Immunosuppressive Niche Score 小提琴图（核心评分，预期在 stroma_immune 域最高） |

---

### Figure 4：签名基因与通路

| 文件名 | 对应图 | 内容 |
|--------|--------|------|
| `fig4A_volcano.png` | Figure 4A | 论文版火山图（niche_high vs niche_low），标注 FOXP3/TGFB1/FAP/CD163/SPP1 等先验基因；高质量输出（300 dpi） |
| `fig4B_niche_marker_heatmap.png` | Figure 4B | Top 30 签名基因在所有 spot 中的表达热图（spot 按 niche_high/niche_low 分组排列），YBP 配色 |
| `fig4C_pathway_bubble.png` | Figure 4C | 签名基因通路富集气泡图（基于 Reactome/KEGG 等通路数据库），气泡大小=基因数，颜色=-log10(p) |

---

### Figure 5：验证与稳健性

| 文件名 | 对应图 | 内容 |
|--------|--------|------|
| `fig5A_label_transfer.png` | Figure 5A | 将 HCC4R 的 niche 签名迁移到 CHC20 切片并可视化（标签迁移验证跨队列可重复性） |
| `fig5B_validation_violin_*.png` | Figure 5B | 跨队列（HCC4R vs CHC20）的各评分小提琴图对比（每个评分一个文件） |
| `fig5C_sensitivity_stability.png` | Figure 5C | 论文版敏感性分析 2×2 图（直接读取 `sensitivity_analysis.csv`） |
| `fig5C_param_scan_heatmap.png` | Figure 5C | 论文版参数扫描热图（直接读取 `param_scan_deg_stability.csv`） |

---

## 五、关键文件依赖关系

```
data/HCC4R/                      ← Visium 原始数据（输入）
data/scRNA_reference.h5ad        ← scRNA 参考数据（输入）
        │
        ▼ Step 1 (run_preprocessing.py)
results/HCC4R/
  ├── adata_vis_post.h5ad         ← 最重要中间结果
  ├── adata_sc_post.h5ad
  ├── spot_cell_proportion_HCC4R.csv
  └── 各种可视化图...
        │
        ▼ Step 2 (run_spatial_niche_analysis.py)
results/HCC4R/spatial_niche/
  ├── spatial_niche_scores.csv    ← 最核心输出
  ├── immunosuppressive_niche_signature_genes_ranked.csv
  ├── immunosuppressive_niche_signature_genes.txt
  └── plots/*.png
results/HCC4R/spatial_signature_genes.txt
        │
        ├──► tcga_survival_analysis.R  (Step 4-5: TCGA ssGSEA 预后分析)
        │
        ▼ Step 3 (run_paper_figures.py，可选)
results/paper_figures/
  └── fig1A ~ fig5C *.png         ← 论文投稿图
```
