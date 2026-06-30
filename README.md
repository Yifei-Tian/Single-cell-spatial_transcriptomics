# 肝癌空间免疫抑制微环境分析流程

本项目基于单细胞转录组（scRNA-seq）与空间转录组（10x Visium）数据，通过 Cell2location 反卷积构建肿瘤免疫抑制空间微环境（Spatial Immunosuppressive Niche）签名，并将其投影到 TCGA-LIHC 队列进行生存预后验证。

**当前分析模式：HCC4R 主分析 + CHC20 独立跨队列验证**
- **HCC4R**（Discovery Cohort）：主分析数据集，运行完整三步流程
- **CHC20**（Validation Cohort）：独立验证数据集，用于 Figure 5 跨队列验证

> **最新更新**：Figure 5 已升级为双模块评分体系（免疫模块 + 基质/ECM 模块），签名基因筛选引入 `signature_rank_score` 综合排序与非特异性基因过滤。详见 [Figure 5 说明](#figure-5-双模块评分验证)。

---

## 目录

- [生物学背景](#生物学背景)
- [分析流程总览](#分析流程总览)
- [代码结构](#代码结构)
- [快速开始](#快速开始)
  - [HCC4R 主分析（推荐）](#hcc4r-主分析推荐)
  - [CHC20 验证分析](#chc20-验证分析)
  - [分步骤运行](#分步骤运行)
  - [TCGA 生存分析（R）](#tcga-生存分析r)
- [输出文件结构](#输出文件结构)
- [Figure 5 双模块评分验证](#figure-5-双模块评分验证)
- [pipeline/ 核心脚本说明](#pipeline-核心脚本说明)
- [utils/ 工具脚本说明](#utils-工具脚本说明)
- [依赖环境](#依赖环境)

---

## 生物学背景

本项目关注的核心生物学问题：

- **肿瘤实质区域**（以 Hepatocyte-like 高丰度 spot 为代表）与**免疫抑制 T 细胞信号**（Treg-like）在空间上是否存在共定位？
- Treg-like 信号是否富集于肿瘤边缘或基质/免疫区域？
- Treg、Myeloid（TAM）、Fibroblast（CAF）三类细胞是否在空间上协同形成**免疫抑制微生态位（Immunosuppressive Spatial Niche）**？
- 是否可以提炼出免疫抑制生态位**基因签名**，并在独立队列（CHC20）和 TCGA-LIHC 大队列中验证其临床预后意义？

---

## 分析流程总览

```
scRNA-seq (GSE149614)             Visium 空间切片 (HCC4R / CHC20)
        |                                        |
   [utils/pre.py / utils/preprocessing.py]       |
   数据预处理 & 质控 & Treg 重注释               质控 & 归一化 & 聚类
        |                                        |
        +-----------> Cell2location <------------+
                  [pipeline/run_preprocessing.py]
                  Step 1：RegressionModel (scRNA) → 参考签名
                         Cell2location (Visium) → spot 细胞类型丰度
                                     |
                  [pipeline/run_spatial_niche_analysis.py]
                  Step 2：kNN 邻域 + Leiden 聚类 → 免疫抑制生态位发现
                          三层基因筛选 → 免疫抑制签名基因
                          参数扫描验证（k=15, q=0.85 为最优参数）
                                     |
                  [pipeline/run_paper_figures.py]（可选）
                  Step 3：生成论文 Figure 1–6 全套高质量子图
                                     |
                  [pipeline/tcga_survival_analysis.R]（可选）
                  Step 4：ssGSEA 投影 → KM + Cox 生存分析
```

---

## 代码结构

```
code/
├── run_hcc4r.py              ← HCC4R 主分析一键入口 ★ 从这里开始
├── run_chc20.py              ← CHC20 独立验证分析入口
│
├── pipeline/                 ← 核心分析流程（被入口脚本调用）
│   ├── run_preprocessing.py       Step 1：数据预处理 + Cell2location 反卷积
│   ├── run_spatial_niche_analysis.py  Step 2：空间免疫抑制生态位分析
│   ├── run_paper_figures.py       Step 3（可选）：论文图表生成（Figure 1–6）
│   ├── paper_plot_functions.py    绘图函数库（PAPER_YBP 配色、DOMAIN_COLORS 等）
│   ├── tcga_survival_analysis.R   Step 4（可选）：TCGA 生存分析
│   └── requirements_r.txt         R 包依赖
│
└── utils/                    ← 工具辅助模块（被 pipeline/ 调用）
    ├── preprocessing.py           核心预处理函数库（Source of Truth）
    ├── pre.py                     早期数据预处理辅助工具（从原始 txt 重建 h5ad）
    ├── colocation.py              逻辑回归共定位评分（可选替代方法）
    └── inspect_data_structure.py  调试辅助工具（查看 AnnData 结构）

data/
├── scRNA_reference.h5ad      scRNA-seq 参考数据（需提前准备）
├── HCC4R/                    HCC4R Space Ranger 输出目录（主分析）
├── HCC1R/                    HCC1R Space Ranger 输出目录（HCC4R 配对验证，可选）
└── CHC20_Visium/             CHC20 Space Ranger 输出目录（独立验证队列）

results/
├── HCC4R/                    HCC4R 主分析结果
│   ├── adata_vis_post.h5ad
│   ├── spatial_niche/
│   └── ...
├── CHC20/                    CHC20 验证分析结果
│   ├── adata_vis_post.h5ad
│   ├── spatial_niche/
│   └── ...
└── paper_figures/            论文图表（Step 3 输出）
    ├── fig1A_workflow_diagram.png
    ├── fig2A_spatial_celltype_all.png
    └── ...（Figure 1–6 全套子图）
```

> **遗留文件说明**：根目录下的 `code/run_preprocessing.py`、`code/preprocessing.py`、`code/run_spatial_niche_analysis.py` 为早期版本，已标记为 Deprecated。请使用 `code/pipeline/` 下的对应脚本。

---

## 快速开始

### HCC4R 主分析（推荐）

```bash
# ① 完整流程：Step 1（反卷积）+ Step 2（生态位分析）
python code/run_hcc4r.py

# ② 完整流程 + 生成论文图表（Figure 1–6）
python code/run_hcc4r.py --paper-figures

# ③ 仅预处理（Step 1）
python code/run_hcc4r.py --step1-only

# ④ 仅空间分析（已有 adata_vis_post.h5ad 时跳过 Step 1）
python code/run_hcc4r.py --step2-only

# ⑤ 仅生成论文图表（Step 1+2 已完成）
python code/run_hcc4r.py --figures-only

# ⑥ Step 1 中跳过 HCC1R 配对验证切片（仅处理 HCC4R 单样本）
python code/run_hcc4r.py --no-sample2
```

**输出目录**：
- `results/HCC4R/` — Cell2location 反卷积 + 生态位分析结果
- `results/paper_figures/` — 论文图表（使用 `--paper-figures` 或 `--figures-only` 时生成）

---

### CHC20 验证分析

CHC20 作为独立验证队列（Validation Cohort），独立运行 Step 1 + Step 2，结果用于论文 Figure 5 跨队列验证。

```bash
# 完整流程（Step 1 + Step 2）
python code/run_chc20.py

# 仅预处理（Step 1）
python code/run_chc20.py --step1-only

# 仅空间分析（Step 2）
python code/run_chc20.py --step2-only
```

**输出目录**：`results/CHC20/`

---

### 分步骤运行

如需直接调用 `pipeline/` 脚本手动控制各步骤：

**Step 1：Cell2location 反卷积**

```bash
# 通过入口脚本调用（推荐）
python code/run_hcc4r.py --step1-only
```

**Step 2：空间生态位分析**

```bash
# HCC4R（主分析，最优参数 k=15, q=0.85 为默认值）
python code/pipeline/run_spatial_niche_analysis.py \
    --adata results/HCC4R/adata_vis_post.h5ad \
    --out-dir results/HCC4R/spatial_niche \
    --signature-out results/HCC4R/spatial_signature_genes.txt

# CHC20（验证分析）
python code/pipeline/run_spatial_niche_analysis.py \
    --adata results/CHC20/adata_vis_post.h5ad \
    --out-dir results/CHC20/spatial_niche \
    --signature-out results/CHC20/spatial_signature_genes.txt
```

`run_spatial_niche_analysis.py` 支持的完整参数：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--adata` | — | 输入 AnnData 路径（必填） |
| `--out-dir` | `results/spatial_niche/` | 输出目录 |
| `--signature-out` | — | 签名基因输出路径 |
| `--n-neighbors` | **15** | kNN 邻居数（参数扫描验证的最优值） |
| `--niche-high-quantile` | **0.85** | niche_high 分位数阈值（参数扫描验证的最优值） |
| `--hep-high-quantile` | `0.75` | 高肝细胞区域分位数阈值 |
| `--leiden-resolution` | `0.5` | Leiden 聚类分辨率 |
| `--top-niche-genes` | `80` | 输出签名基因数量 |
| `--per-sample-neighbors` | False | 联合分析模式：kNN 邻域只在同一切片内建立 |

**Step 3：论文图表生成**

```bash
# 生成全套论文图表（Figure 1–6）
python code/run_hcc4r.py --figures-only

# 或直接调用
python code/pipeline/run_paper_figures.py \
    --hcc4r-dir results/HCC4R \
    --chc20-dir results/CHC20 \
    --out-dir results/paper_figures
```

> **修改绘图样式**：编辑 `code/pipeline/paper_plot_functions.py`，无需改动数据准备逻辑。

---

### TCGA 生存分析（R）

基于 Step 2 生成的签名基因列表，在 TCGA-LIHC 队列中进行 ssGSEA 评分 + KM + Cox 生存分析：

```bash
# 安装 R 依赖
Rscript -e "install.packages(c('GSVA','survival','survminer','ggplot2','forestplot'))"

# 使用 HCC4R 的签名基因（推荐）
Rscript code/pipeline/tcga_survival_analysis.R \
    --signature results/HCC4R/spatial_niche/immunosuppressive_niche_signature_genes.txt \
    --out-dir results/HCC4R/tcga

# 使用 CHC20 的签名基因
Rscript code/pipeline/tcga_survival_analysis.R \
    --signature results/CHC20/spatial_niche/immunosuppressive_niche_signature_genes.txt \
    --out-dir results/CHC20/tcga
```

---

## 输出文件结构

每个数据集的完整输出目录（以 `results/HCC4R/` 为例）：

```
results/HCC4R/
├── adata_vis_post.h5ad                        - Cell2location 反卷积后的 Visium AnnData（主文件）
├── adata_sc_post.h5ad                         - RegressionModel 训练后的 scRNA AnnData
├── spot_cell_proportion_HCC4R.csv             - 每个 spot 的细胞类型归一化比例
├── shared_genes_HCC4R.txt                     - scRNA ∩ Visium 共有基因列表
├── regression_training_history_HCC4R.png      - RegressionModel 训练损失曲线
├── spatial_mapping_training_history_HCC4R.png - Cell2location 空间建模训练曲线
├── scrna_tsne_celltype.png                    - scRNA tSNE 细胞类型图（论文 Fig 1B）
├── scrna_celltype_marker_heatmap.png          - Marker 基因热图（论文 Fig 1C）
├── t_cell_dotplot_horizontal.png              - Treg 标志基因气泡图（论文 Fig 1D）
├── spatial_signature_genes.txt                - 最终签名基因（TCGA 分析接口）
│
├── spatial_niche/                             - Step 2 空间生态位分析结果
│   ├── spatial_niche_scores.csv              - 每个 spot 的全部评分指标（核心数据表）
│   ├── spatial_niche_parameters.csv          - 分析参数记录（k=15, q=0.85 等）
│   ├── neighborhood_cluster_stats.csv        - Leiden 聚类各维度统计
│   ├── sensitivity_analysis.csv              - kNN/radius 敏感性分析（连续指标）
│   ├── immunosuppressive_niche_signature_genes_ranked.csv  - 三层筛选签名基因排名表
│   ├── immunosuppressive_niche_signature_genes.txt         - 签名基因列表（TCGA 接口）
│   ├── prior_gene_set_auc.csv               - Layer 2 先验基因集 AUC 结果
│   ├── gini_score_genes.csv                 - Layer 3 Gini 特异性基因（阈值 Gini>0.3）
│   ├── param_scan_deg_stability.csv         - k × quantile 参数扫描结果（验证 k=15, q=0.85）
│   └── plots/                               - 全套空间可视化图（paper_ybp 配色）
│       ├── spatial_hepatocyte.png
│       ├── spatial_treg.png
│       ├── spatial_myeloid.png
│       ├── spatial_fibroblast.png
│       ├── spatial_immunosuppressive_niche_score.png
│       ├── spatial_neighborhood_clusters.png
│       ├── spatial_niche_semantic_labels.png
│       ├── spatial_region_labels.png
│       ├── spatial_niche_high_score_spots.png
│       ├── spatial_niche_cluster_vs_score_comparison.png
│       ├── niche_signature_volcano.png
│       ├── niche_fraction_scatter.png
│       ├── lr_communication_heatmap.png
│       ├── celltype_niche_correlation.png
│       ├── hepatocyte_vs_treg_niche_score.png
│       ├── distance_to_hep_high_vs_niche_score.png
│       ├── region_score_boxplots.png
│       ├── sensitivity_niche_stability.png
│       └── param_scan_deg_stability_heatmap.png
│
└── cross_slice_comparison/                    - 配对验证切片（HCC1R）一致性对比图
    ├── cross_slice_mean_proportion_comparison.png
    ├── cross_slice_treg_distribution.png
    └── cross_slice_celltype_boxplot.png

results/paper_figures/                         - Step 3 论文图表（DPI=300，期刊投稿质量）
├── fig1A_workflow_diagram.png                 - 课题技术路线图
├── fig1B_scrna_tsne_celltype.png              - scRNA tSNE 细胞类型图
├── fig1C_scrna_marker_heatmap.png             - Marker 基因热图
├── fig1D_t_cell_dotplot.png                   - Treg 标志基因气泡图
├── fig2A_spatial_celltype_all.png             - 所有细胞类型 3×3 空间分布网格
├── fig2B_spatial_hepatocyte.png               - Hepatocyte 空间分布
├── fig2C_spatial_treg.png                     - Treg 空间分布
├── fig2D_spatial_myeloid.png                  - Myeloid 空间分布
├── fig2E_spatial_fibroblast.png               - Fibroblast 空间分布
├── fig3A_spatial_niche_clusters.png           - Leiden 邻域聚类空间图
├── fig3B_spatial_niche_semantic.png           - 语义 niche 标签图
├── fig3C_spatial_niche_score.png              - 综合 niche 评分空间热图
├── fig3D_spatial_niche_high.png               - niche_high（q=0.85）二值分布图
├── fig3E_region_score_boxplots.png            - 各区域评分箱线图
├── fig3F_hepatocyte_treg_scatter.png          - Hepatocyte vs Treg 散点图
├── fig4A_niche_volcano.png                    - 签名基因火山图
├── fig4B_niche_fraction_scatter.png           - 检出率差值散点图（捕获 FOXP3 等稀有基因）
├── fig4C_lr_communication_heatmap.png         - 配体-受体通讯分析热图
├── fig4D_celltype_correlation.png             - 细胞类型相关性热图
├── fig5A_chc20_immune_module_score.png        - CHC20 免疫模块空间评分图（Fig.5A 新版）
├── fig5B_chc20_stromal_ECM_module_score.png   - CHC20 基质/ECM 模块空间评分图（Fig.5B 新版）
├── fig5C_chc20_combined_niche_score.png       - CHC20 组合 niche score（z_immune + z_stromal，Fig.5C）
├── fig5D_combined_score_HCC4R_vs_CHC20.png   - HCC4R vs CHC20 combined niche score 分布对比（Fig.5D）
├── fig5_dual_module_scores_HCC4R.csv          - HCC4R 三列评分数据表（免疫/基质/组合）
├── fig5_dual_module_scores_CHC20.csv          - CHC20 三列评分数据表（免疫/基质/组合）
├── fig5_immune_genes_used_CHC20.txt           - CHC20 实际匹配到的免疫模块基因列表
├── fig5_stromal_genes_used_CHC20.txt          - CHC20 实际匹配到的基质模块基因列表
└── fig5C_sensitivity_stability.png            - 敏感性分析（Spearman ρ + Jaccard，补充图）
```

> 输出文件详细说明见 [`docs/result_explain.md`](docs/result_explain.md)

---

## Figure 5 双模块评分验证

Figure 5 对应论文中的 **CHC20 跨队列验证部分**，已从单一签名评分升级为双模块评分体系，更准确地反映免疫抑制空间微环境的双重特征（免疫抑制信号 + 基质重塑信号）。

### 评分体系

| 评分名称 | 计算公式 | 生物学含义 |
|---|---|---|
| `immune_signature_score` | mean(E_ig, g ∈ 免疫模块基因) | 反映 Treg/TAM 介导的免疫抑制强度 |
| `stromal_ECM_signature_score` | mean(E_ig, g ∈ 基质模块基因) | 反映 CAF/ECM 重塑介导的屏障效应 |
| `immune_stromal_niche_score` | z(immune) + z(stromal) | 综合 niche 强度（两模块 Z-score 之和） |

其中 `E_ig` 为 log1p-normalized(10k) 表达量，`z(·)` 为跨 spot 的 Z-score 标准化。

### 生成的图表（输出到 `results/paper_figures/`）

| 图表文件 | 内容 |
|---|---|
| `fig5A_chc20_immune_module_score.png` | CHC20 免疫模块空间投影（paper_ybp 配色） |
| `fig5B_chc20_stromal_ECM_module_score.png` | CHC20 基质/ECM 模块空间投影 |
| `fig5C_chc20_combined_niche_score.png` | CHC20 组合 niche score 空间投影 |
| `fig5D_combined_score_HCC4R_vs_CHC20.png` | HCC4R vs CHC20 评分分布对比（Mann-Whitney U 检验） |

### 自定义模块基因

**方式 1**（推荐）：在 `immunosuppressive_niche_signature_genes_ranked.csv` 中增加 `module` 列，取值为 `"immune"` 或 `"stromal"`/`"ecm"`，脚本将自动读取。

**方式 2**：直接修改 `run_paper_figures.py` 中 `_split_module_genes()` 函数内的 `_IMMUNE_GENES_DEFAULT` 和 `_STROMAL_GENES_DEFAULT` 列表。

---

## pipeline/ 核心脚本说明

### `pipeline/run_preprocessing.py`

**功能**：Step 1 主流程，对一张（或两张）Visium 切片完成 Cell2location 反卷积。

| 内部步骤 | 说明 |
|---|---|
| Step 1–3 | 加载 scRNA / Visium 数据，取共同基因子集 |
| Step 4 | 提取 T/NK 亚群，将簇 9 重注释为 Treg（FOXP3+CD4+） |
| Step 5 | 训练 RegressionModel（学习各细胞类型参考签名） |
| Step 6 | 运行 Cell2location 空间建模（主切片） |
| Step 7–8 | 验证切片独立建模 + 多切片一致性对比图（有 sample2 时执行） |
| Step 9 | 保存共享基因列表和所有 AnnData / 可视化图 |

`main()` 函数接受 `path_sample1`, `sample1_name`, `output_dir` 等参数，不同数据集只需传入不同参数，无需修改脚本内容。

---

### `pipeline/run_spatial_niche_analysis.py`

**功能**：Step 2 主流程，基于 Cell2location 反卷积结果定义免疫抑制空间生态位，提取生态位特异性基因签名。

**两阶段核心逻辑**：

1. **无监督发现**：k=15 kNN 邻域平滑 → Leiden 聚类 → 四维评分排名 → 识别 `immunosuppressive_niche`
2. **基因筛选（三层，已升级）**：
   - Layer 1：Wilcoxon + BH-FDR + AUC + `signature_rank_score` 综合排序
     - 筛选条件：`FDR < 0.1 AND AUC > 0.60 AND (log2FC > 0.3 OR delta_frac > 0.10)`
     - 自动放宽：基因数 < 80 时，放宽至 `FDR < 0.2 AND AUC > 0.55 AND (FC > 0.2 OR Δfrac > 0.05)`
     - 非特异性基因过滤：管家基因、核糖体基因、线粒体基因、血红蛋白基因、细胞周期基因等
   - Layer 2：先验基因集 AUC（Treg/TAM/CAF 先验基因强制纳入，排除非特异性基因）
   - Layer 3：Gini Index（阈值 0.3，捕获 FOXP3 等稀有局灶性基因）

**签名基因综合排序公式**：

```
signature_rank_score = 0.35 × minmax(log2FC)
                     + 0.35 × minmax(delta_frac)
                     + 0.30 × minmax(AUC)
```

**关键设计**：

| 特性 | 说明 |
|---|---|
| 最优参数 k=15, q=0.85 | 经 k × quantile 二维参数扫描验证（`param_scan_deg_stability_heatmap.png`） |
| delta_frac 指标 | 计算 spot 检出率差值（`frac_high - frac_low`），解决 Visium 稀释效应对稀有基因（FOXP3）不敏感的问题 |
| signature_rank_score | 三维综合排序评分（log2FC + delta_frac + AUC 加权融合） |
| 非特异性基因过滤 | 管家/核糖体/线粒体/血红蛋白/细胞周期/广谱炎症/组织损伤基因黑名单（共 350+ 个） |
| Gini Index 阈值 0.3 | 覆盖局灶性高表达稀有免疫基因（原为 0.5，已优化） |
| `--per-sample-neighbors` | 联合分析时可启用，kNN 邻域只在同一切片内建立 |

---

### `pipeline/run_paper_figures.py`

**功能**：Step 3，生成适合期刊投稿的高质量论文图表（DPI=300）。

- 基于 Step 1（HCC4R AnnData）和 Step 2（niche scores CSV）的输出
- 生成论文 Figure 1–5 的所有子图，每张独立保存到 `results/paper_figures/`
- Figure 5 需要 CHC20 的分析结果（`results/CHC20/adata_vis_post.h5ad`）

**Figure 5 新版流程（双模块评分）**：

1. 从 `immunosuppressive_niche_signature_genes_ranked.csv` 读取 `module` 列拆分免疫/基质基因
2. 若无 `module` 列，自动使用内置先验基因列表（免疫模块 21 个 + 基质模块 21 个）
3. 分别在 HCC4R 和 CHC20 上计算三套评分：
   - `immune_signature_score` = 免疫模块基因均值表达
   - `stromal_ECM_signature_score` = 基质/ECM 模块基因均值表达
   - `immune_stromal_niche_score` = z(immune) + z(stromal)（跨 spot Z-score 之和）
4. 输出 Fig.5A（免疫模块空间图）、Fig.5B（基质模块空间图）、Fig.5C（组合 niche 空间图）、Fig.5D（HCC4R vs CHC20 对比小提琴图）

**修改绘图样式**：只需编辑 `code/pipeline/paper_plot_functions.py`（绘图函数库），无需改动数据准备逻辑：

```python
# paper_plot_functions.py 中的关键常量
PAPER_YBP = ...          # Yellow-Black-Purple 三色渐变色标（空间热图主配色）
DOMAIN_COLORS = {...}    # 各细胞类型/区域的标准颜色映射
```

---

### `pipeline/tcga_survival_analysis.R`

**功能**：Step 4（可选），将签名基因集投影到 TCGA-LIHC 大队列（n≈370），进行 ssGSEA 评分 + Kaplan-Meier + Cox 多变量生存分析。

**主要输出**：

| 文件 | 内容 |
|---|---|
| `tcga_signature_score.csv` | 每位患者的 ssGSEA 免疫抑制评分 |
| `tcga_signature_survival.csv` | 评分 + 临床信息合并表 |
| `cox_results.txt` | Cox 比例风险回归结果 |
| `km_plot.png` | Kaplan-Meier 生存曲线（高/低评分组） |
| `cox_forest_plot.png` | 多变量 Cox 森林图 |

---

## utils/ 工具脚本说明

### `utils/preprocessing.py`

核心预处理函数库（Source of Truth），被 `pipeline/run_preprocessing.py` 调用。封装了从 h5ad 加载数据到 Cell2location 建模的所有中间步骤（`load_scrna_h5ad`, `load_visium`, `align_shared_genes`, `assign_treg_label`, `setup_and_train_regression_model` 等）。

### `utils/pre.py`

早期数据预处理辅助脚本，用于从 txt 格式原始数据重建 scRNA-seq 参考 h5ad。如需从头重建 `scRNA_reference.h5ad`：

```bash
python code/utils/pre.py
```

### `utils/colocation.py`

独立的共定位辅助分析脚本，提供基于逻辑回归的软性共定位评分方法（可选替代方案，供探索性分析使用）。

### `utils/inspect_data_structure.py`

调试辅助工具，快速检查 AnnData 对象的数据结构：

```bash
python code/utils/inspect_data_structure.py \
    --sc-h5ad results/HCC4R/adata_sc_post.h5ad \
    --spatial-h5ad results/HCC4R/adata_vis_post.h5ad \
    --query sc_obs --head 10
```

---

## 依赖环境

### Python 依赖

```bash
pip install -r requirements.txt
```

主要依赖：

```
scanpy
anndata
squidpy
cell2location
scipy
numpy
pandas
seaborn
matplotlib
scikit-learn
torch
statsmodels
```

### R 依赖

```r
install.packages(c("survival", "survminer", "ggplot2", "forestplot", "dplyr"))
BiocManager::install("GSVA")
```

详见 `code/pipeline/requirements_r.txt`
