# 肝癌空间免疫抑制微环境分析流程

本项目基于单细胞转录组（scRNA-seq）与空间转录组（10x Visium）数据，通过 Cell2location 反卷积，构建肿瘤免疫抑制空间微环境（Spatial Immunosuppressive Niche）签名，并将其投影到 TCGA-LIHC 队列进行生存预后验证。

---

## 目录

- [生物学背景](#生物学背景)
- [分析流程总览](#分析流程总览)
- [代码结构](#代码结构)
- [各数据集运行教程](#各数据集运行教程)
  - [快速开始（推荐）](#快速开始推荐)
  - [CHC20 数据集](#chc20-数据集)
  - [HCC4R 数据集](#hcc4r-数据集)
  - [HCC6NR 数据集](#hcc6nr-数据集)
  - [HCC4R+CHC20 联合分析](#hcc4rchc20-联合分析)
  - [分步骤运行](#分步骤运行)
  - [Step 3（可选）：Wilcoxon 签名基因交叉验证](#step-3可选wilcoxon-签名基因交叉验证)
  - [Step 4 & 5：TCGA 生存分析（R）](#step-4--5tcga-生存分析r)
- [输出文件结构](#输出文件结构)
- [pipeline/ 核心脚本说明](#pipeline-核心脚本说明)
- [utils/ 工具脚本说明](#utils-工具脚本说明)
- [当前代码中的不合理之处与改进建议](#当前代码中的不合理之处与改进建议)

---

## 生物学背景

本项目关注的核心生物学问题：

- **肿瘤实质区域**（以 Hepatocyte-like 高丰度 spot 为代表）与**免疫抑制 T 细胞信号**（Treg-like）在空间上是否存在共定位？
- Treg-like 信号是否富集于肿瘤边缘或基质/免疫区域？
- T/NK、Myeloid、Fibroblast 与 Hepatocyte-like 区域在空间上如何分布和相互关系？
- 是否可以提炼出一个**免疫抑制空间生态位（Immunosuppressive Spatial Niche）基因签名**，并在 TCGA-LIHC bulk RNA-seq 中验证其预后意义？

---

## 分析流程总览

```
scRNA-seq (GSE149614)          Visium 空间切片 (CHC20 / HCC4R / HCC6NR)
        |                                  |
   [utils/pre.py / utils/preprocessing.py] |
   数据预处理 & 质控                        质控 & 归一化 & 聚类
        |                                  |
        +---------> Cell2location <---------+
                  [pipeline/run_preprocessing.py]
                  RegressionModel (scRNA) → 细胞类型参考签名
                  Cell2location (Visium) → 每个 spot 的细胞类型丰度
                         |
              [pipeline/run_spatial_niche_analysis.py]
              构建免疫抑制空间生态位评分 & 基因签名
                         |
              [pipeline/run_de_analysis.py] (可选 Wilcoxon 交叉验证)
                         |
              [pipeline/tcga_survival_analysis.R]
              ssGSEA 投影 → KM + Cox 生存分析
```

---

## 代码结构

```
code/
├── run_chc20.py              ← CHC20 数据集一键入口（直接运行这个）
├── run_hcc4r.py              ← HCC4R 数据集一键入口（直接运行这个）
├── run_hcc6nr.py             ← HCC6NR 数据集一键入口（直接运行这个）
├── run_joint_hcc4r_chc20.py  ← HCC4R+CHC20 联合分析入口（直接运行这个）
│
├── pipeline/                 ← 核心分析流程（被入口脚本调用）
│   ├── run_preprocessing.py      Step 1：数据预处理 + Cell2location 反卷积
│   ├── run_spatial_niche_analysis.py  Step 2：空间免疫抑制生态位分析
│   ├── run_de_analysis.py        Step 3（可选）：Wilcoxon 差异表达验证
│   ├── run_chc23_validation.py   CHC23 切片独立验证（历史脚本）
│   ├── tcga_survival_analysis.R  Step 4-5：TCGA 生存分析
│   └── requirements_r.txt        R 包依赖
│
└── utils/                    ← 工具辅助模块（被 pipeline/ 调用）
    ├── preprocessing.py          核心预处理函数库
    ├── pre.py                    早期数据预处理辅助工具
    ├── colocation.py             逻辑回归共定位评分（可选替代方法）
    └── inspect_data_structure.py 调试辅助工具（查看 AnnData 结构）

data/
├── scRNA_reference.h5ad      scRNA-seq 参考数据（需提前准备）
├── CHC20_Visium/             CHC20 Space Ranger 输出目录
├── CHC23_Visium/             CHC23 Space Ranger 输出目录
├── HCC4R/                    HCC4R Space Ranger 输出目录（新增）
└── HCC6NR/                   HCC6NR Space Ranger 输出目录（新增）

results/
├── CHC20/                    CHC20 所有分析结果
│   ├── adata_vis_post.h5ad
│   ├── spatial_niche/
│   └── ...
├── HCC4R/                    HCC4R 所有分析结果
│   ├── adata_vis_post.h5ad
│   ├── spatial_niche/
│   └── ...
├── HCC6NR/                   HCC6NR 所有分析结果
│   ├── adata_vis_post.h5ad
│   ├── spatial_niche/
│   └── ...
└── joint_HCC4R_CHC20/        HCC4R+CHC20 联合分析结果
    ├── adata_vis_post.h5ad           合并反卷积结果（含 obs["sample"] 列）
    ├── adata_vis_post_HCC4R.h5ad     HCC4R 单独结果
    ├── adata_vis_post_CHC20.h5ad     CHC20 单独结果
    ├── spot_cell_proportion_joint.csv 合并比例表
    ├── spatial_signature_genes.txt   联合分析签名基因
    ├── spatial_niche/                联合空间生态位分析结果
    │   └── ...（同单数据集结构）
    └── cross_slice_comparison/       两切片一致性对比图
```

---

## 各数据集运行教程

### 快速开始（推荐）

每个数据集只需运行对应的入口脚本，**无需修改任何代码**：

```bash
# CHC20（原始主分析数据集）
python code/run_chc20.py

# HCC4R（新数据集，验证效果）
python code/run_hcc4r.py

# HCC6NR（新数据集，验证效果）
python code/run_hcc6nr.py
```

每个入口脚本会自动：
1. 执行 Step 1（预处理 + Cell2location 反卷积）
2. 执行 Step 2（空间生态位分析）
3. 将结果保存到各自的 `results/<数据集名>/` 目录

如需跨数据集联合分析（如 HCC4R+CHC20），请使用 [联合分析入口](#hcc4rchc20-联合分析)。

---

### CHC20 数据集

```bash
# 完整流程（Step 1 + Step 2）
python code/run_chc20.py

# 仅预处理（Step 1）
python code/run_chc20.py --step1-only

# 仅空间分析（Step 2，已有 adata_vis_post.h5ad 时使用）
python code/run_chc20.py --step2-only
```

**输出目录**：`results/CHC20/`

---

### HCC4R 数据集

```bash
# 完整流程（Step 1 + Step 2）
python code/run_hcc4r.py

# 仅预处理
python code/run_hcc4r.py --step1-only

# 仅空间分析
python code/run_hcc4r.py --step2-only

# Step 1 中跳过 HCC6NR 配对验证切片（仅处理 HCC4R 单样本）
python code/run_hcc4r.py --no-sample2
```

**输出目录**：`results/HCC4R/`

---

### HCC6NR 数据集

```bash
# 完整流程（Step 1 + Step 2）
python code/run_hcc6nr.py

# 仅预处理
python code/run_hcc6nr.py --step1-only

# 仅空间分析
python code/run_hcc6nr.py --step2-only

# Step 1 中跳过 HCC4R 配对验证切片（仅处理 HCC6NR 单样本）
python code/run_hcc6nr.py --no-sample2
```

**输出目录**：`results/HCC6NR/`

---

### HCC4R+CHC20 联合分析

经批次效应评估确认 HCC4R 与 CHC20 的 batch effect 极小，适合联合分析。联合分析使用同一套 scRNA 参考签名训练一次 RegressionModel，分别反卷积后合并 spot 级结果，在合并数据上进行空间生态位识别（邻域计算限制在切片内）。

```bash
# 完整流程（Step 1 + Step 2）
python code/run_joint_hcc4r_chc20.py

# 仅预处理 + 反卷积 + 合并
python code/run_joint_hcc4r_chc20.py --step1-only

# 仅空间生态位联合分析（已有 adata_vis_post.h5ad 时）
python code/run_joint_hcc4r_chc20.py --step2-only
```

**关键参数**：`--per-sample-neighbors`（Step 2 自动启用）
- 空间 kNN 邻域只在同一切片内建立，不允许 HCC4R 和 CHC20 的 spot 跨切片互为物理邻居
- 需要合并后的 `adata.obs['sample']` 列标记切片来源
- 下游 cell composition、cluster、niche score、DEG 分析在合并后的全表上进行，获得更大样本量和更高统计功效

**输出目录**：`results/joint_HCC4R_CHC20/`

---

### 分步骤运行

如需对某个数据集手动控制各步骤，也可以直接调用 `pipeline/` 中的脚本：

**Step 1：Cell2location 反卷积**

```bash
# 直接运行（使用默认 CHC20/CHC23 路径）
python code/pipeline/run_preprocessing.py
```

**Step 2：空间生态位分析**

```bash
# CHC20
python code/pipeline/run_spatial_niche_analysis.py \
    --adata results/CHC20/adata_vis_post.h5ad \
    --out-dir results/CHC20/spatial_niche \
    --signature-out results/CHC20/spatial_signature_genes.txt

# HCC4R
python code/pipeline/run_spatial_niche_analysis.py \
    --adata results/HCC4R/adata_vis_post.h5ad \
    --out-dir results/HCC4R/spatial_niche \
    --signature-out results/HCC4R/spatial_signature_genes.txt

# HCC6NR
python code/pipeline/run_spatial_niche_analysis.py \
    --adata results/HCC6NR/adata_vis_post.h5ad \
    --out-dir results/HCC6NR/spatial_niche \
    --signature-out results/HCC6NR/spatial_signature_genes.txt
```

支持的可选参数：

```bash
python code/pipeline/run_spatial_niche_analysis.py \
    --adata results/HCC4R/adata_vis_post.h5ad \
    --out-dir results/HCC4R/spatial_niche \
    --hep-high-quantile 0.75 \
    --niche-high-quantile 0.80 \
    --top-niche-genes 80 \
    --per-sample-neighbors   # 联合分析时启用：限制邻域在切片内建立
```

| 参数 | 说明 |
|------|------|
| `--per-sample-neighbors` | 联合分析模式下，kNN 邻域只在同一切片（sample）内建立，不允许跨切片互为物理邻居。需要 `adata.obs['sample']` 列 |
| `--hep-high-quantile` | 高肝细胞区域分位数阈值（默认 0.75） |
| `--niche-high-quantile` | 免疫抑制 niche 高分组分位数阈值（默认 0.80） |
| `--top-niche-genes` | 输出 niche 特征基因数量（默认 80） |
| `--n-neighbors` | kNN 邻居数（默认 15） |
| `--leiden-resolution` | Leiden 聚类分辨率（默认 0.5） |

---

### Step 3（可选）：Wilcoxon 签名基因交叉验证

```bash
# HCC4R 示例
python code/pipeline/run_de_analysis.py \
    --coloc results/HCC4R/spatial_niche/spatial_niche_scores.csv \
    --coloc-column niche_high \
    --out results/HCC4R/spatial_niche/wilcoxon_niche_signature_genes.txt \
    --top-n 80
```

---

### Step 4 & 5：TCGA 生存分析（R）

```bash
# 使用 HCC4R 的签名基因
Rscript code/pipeline/tcga_survival_analysis.R \
    --signature results/HCC4R/spatial_signature_genes.txt \
    --out-dir results/HCC4R/tcga

# 使用 CHC20 的签名基因
Rscript code/pipeline/tcga_survival_analysis.R \
    --signature results/CHC20/spatial_signature_genes.txt \
    --out-dir results/CHC20/tcga
```

---

## 输出文件结构

每个数据集的输出目录结构（以 `results/HCC4R/` 为例）：

```
results/HCC4R/
├── run_preprocessing.log                 - 预处理全流程日志
├── adata_vis_post.h5ad                   - Cell2location 反卷积后的 Visium AnnData（主）
├── adata_vis_post_HCC4R.h5ad             - 同上（带数据集名后缀）
├── adata_sc_post.h5ad                    - 训练后的 scRNA AnnData（主）
├── adata_sc_post_HCC4R.h5ad             - 同上（带后缀）
├── spot_cell_proportion_HCC4R.csv        - 每个 spot 的细胞类型比例
├── shared_genes_HCC4R.txt                - scRNA × Visium 共享基因列表
├── regression_training_history_HCC4R.png - RegressionModel 训练曲线
├── spatial_mapping_training_history_HCC4R.png - Cell2location 训练曲线
├── t_cell_dotplot_horizontal.png         - Treg 标志基因气泡图
├── spatial_signature_genes.txt           - 最终签名基因（供 TCGA 分析）
│
├── spatial_niche/                        - Step 2 空间生态位分析结果
│   ├── spatial_niche_scores.csv          - 每个 spot 的生态位评分
│   ├── spatial_niche_parameters.csv      - 分析参数记录
│   ├── immunosuppressive_niche_signature_genes.txt
│   ├── immunosuppressive_niche_signature_genes_ranked.csv
│   ├── gini_score_genes.csv              - Gini Index 特异性基因（Gini>0.3）
│   ├── deg_results_wilcoxon.csv          - DEG 统计结果（含 delta_frac 列）
│   ├── param_scan_deg_stability.csv      - k × quantile 参数扫描结果
│   └── plots/                            - 可视化图表
│       ├── spatial_hepatocyte.png
│       ├── spatial_niche_score.png
│       ├── niche_signature_volcano.png
│       ├── niche_fraction_scatter.png
│       ├── param_scan_deg_stability_heatmap.png
│       └── ...（共约 10 张图）
│
└── cross_slice_comparison/               - 多切片一致性对比图（有配对切片时生成）
    ├── cross_slice_mean_proportion_comparison.png
    ├── cross_slice_treg_distribution.png
    └── cross_slice_celltype_boxplot.png
```

---

## pipeline/ 核心脚本说明

### `pipeline/run_preprocessing.py`

**功能定位**：整个分析的 **Step 1 主流程**，对一张（或两张）Visium 切片完成 Cell2location 反卷积。

**主要步骤**：

| 步骤 | 说明 |
|---|---|
| Step 1 | 加载 scRNA-seq 参考数据（质控过滤） |
| Step 2 | 加载 Visium 空间数据（spot 质控、归一化、聚类） |
| Step 3 | scRNA 与 Visium 基因对齐（取共同基因子集） |
| Step 4 | 提取 T/NK 亚群，将簇 9 重注释为 Treg |
| Step 5 | 训练 RegressionModel（主切片参考签名） |
| Step 6 | 运行 Cell2location 空间建模（主切片） |
| Step 7 | 验证切片独立建模（可选，有 sample2 时执行） |
| Step 8 | 多切片一致性对比图（有 sample2 时绘制） |
| Step 9 | 保存共享基因列表和所有 AnnData |

**参数化设计（新）**：`main()` 函数接受 `path_sample1`, `sample1_name`, `output_dir` 等参数，
不同数据集只需传入不同参数，无需修改脚本内容。直接运行时默认使用 CHC20/CHC23。

**联合分析模式（新）**：`main_joint()` 函数支持两切片联合反卷积：
- 用同一个 scRNA 参考训练一次 RegressionModel，得到统一 `cell_state_df`
- 分别对两张 Visium 切片做 Cell2location 反卷积
- 合并 spot 级结果，添加 `obs["sample"]` 列标记来源
- 入口：`python code/run_joint_hcc4r_chc20.py --step1-only`

---

### `pipeline/run_spatial_niche_analysis.py`

**功能定位**：整个分析的 **Step 2 主运行脚本**，基于 Cell2location 反卷积结果定义免疫抑制空间生态位评分，提取生态位特异性基因签名。

**关键改进（小院士修改）**：

| 改进项 | 详情 |
|---|---|
| DEG 稳定性参数扫描（Step 15） | 从 resolution × quantile 改为 **k × quantile** 二维扫描，k 直接影响 niche_score，热图有实际意义 |
| Gini Index 阈值 | 从 0.5 降至 **0.3**，覆盖 FOXP3 等稀有基因 |
| gini_score_genes.csv | **始终写出**（即使结果为空），保证文件路径可预期 |
| delta_frac 检出率 | 计算 niche_high vs niche_low 的 spot 检出率差值，量化稀疏基因富集 |
| `--per-sample-neighbors` | 联合分析模式下，kNN 邻域只在同一切片内建立，不允许跨切片互为物理邻居 |
| `_build_knn_neighbors_per_sample()` | 按切片分组构建 kNN 邻域，确保空间邻域的生物学意义 |
| `_spatial_neighbors_per_sample()` | 按切片分组构建半径邻域，用于邻域均值特征计算 |

**主要输出**：

| 文件 | 说明 |
|---|---|
| `spatial_niche_scores.csv` | 每个 spot 的全部评分指标 |
| `immunosuppressive_niche_signature_genes_ranked.csv` | 带 log2FC / delta_frac 的完整签名 |
| `gini_score_genes.csv` | Gini Index 特异性基因（Gini>0.3, log2FC>0） |
| `param_scan_deg_stability.csv` | k × quantile 参数扫描结果 |
| `plots/param_scan_deg_stability_heatmap.png` | 参数扫描热图（颜色深=DEG 数量多且稳定） |

---

### `pipeline/run_de_analysis.py`

**功能定位**：**可选的 Step 3**，使用 Scanpy Wilcoxon 秩和检验对 niche_high/niche_low 两组 spot 做差异表达验证，提供签名基因的独立交叉验证。

---

### `pipeline/tcga_survival_analysis.R`

**功能定位**：**Step 4 & 5**，将空间生态位签名投影到 TCGA-LIHC 队列，进行 ssGSEA 评分 + KM + Cox 生存分析。

---

## utils/ 工具脚本说明

### `utils/preprocessing.py`

核心预处理函数库，被 `pipeline/run_preprocessing.py` 调用。封装了从 h5ad 加载数据到 Cell2location 建模的所有中间步骤（`load_scrna_h5ad`, `load_visium`, `align_shared_genes`, `assign_treg_label`, `setup_and_train_regression_model` 等）。

### `utils/pre.py`

早期数据预处理辅助脚本（探索性工具），用于从 txt 格式构建 scRNA-seq 参考数据，及加载合并多 Visium 切片。如需从原始数据重建 `scRNA_reference.h5ad`，运行：

```bash
python code/utils/pre.py
```

### `utils/colocation.py`

独立的共定位辅助分析脚本，提供基于逻辑回归的软性共定位评分方法（可选替代方案）。

### `utils/inspect_data_structure.py`

调试辅助工具，快速检查 AnnData 对象的数据结构：

```bash
python code/utils/inspect_data_structure.py \
    --sc-h5ad results/HCC4R/adata_sc_post.h5ad \
    --spatial-h5ad results/HCC4R/adata_vis_post.h5ad \
    --query sc_obs --head 10
```

---

## 当前代码中的不合理之处与改进建议

### 1. 两套预处理模块并存（高优先级）

**问题**：`utils/pre.py` 和 `utils/preprocessing.py` 存在大量功能重叠，且基因对齐方式不一致（前者用 Ensembl ID，后者用 Gene Symbol）。

**建议**：合并为单一预处理模块，明确以 Ensembl ID 为统一基因索引标准。

---

### 2. `q05` 与 `means` 丰度选择缺乏明确说明（中优先级）

**问题**：`run_preprocessing.py` 历史版本注释中有 `# TODO: 这里太保守`，而 `run_spatial_niche_analysis.py` 默认使用 `means_cell_abundance_w_sf`。

**建议**：统一使用 `means_cell_abundance_w_sf` 作为默认丰度指标，q05 保留为保守估计备用选项。

---

### 3. RegressionModel 验证切片训练轮次设置（中优先级）

**当前设置**：验证切片（CHC23/HCC6NR 等）使用 `max_epochs=400, patience=50, min_delta=5e-5`，比主切片更保守，以防早停过早。若新数据集表现不佳，可进一步调整。

---

### 4. `colocation.py` 的逻辑回归方法存在数据泄露风险（低优先级）

**问题**：训练目标直接由特征线性组合而来，模型学到的本质上是同一线性边界，不具真正泛化价值。

**建议**：改用无监督方法（GMM、DBSCAN）或加入更多生物学特征。

---

### 5. 缺少统一的 Python 环境依赖文件（低优先级）

**建议**：参考项目根目录的 `requirements.txt`，确保所有依赖已列出：

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