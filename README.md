# 肝癌空间免疫抑制微环境分析流程

本项目基于单细胞转录组（scRNA-seq）与空间转录组（10x Visium）数据，通过 Cell2location 反卷积，构建肿瘤免疫抑制空间微环境（Spatial Immunosuppressive Niche）签名，并将其投影到 TCGA-LIHC 队列进行生存预后验证。

---

## 目录

- [生物学背景](#生物学背景)
- [分析流程总览](#分析流程总览)
- [文件结构说明](#文件结构说明)
  - [pre.py](#prepy)
  - [preprocessing.py](#preprocessingpy)
  - [run_preprocessing.py](#run_preprocessingpy)
  - [run_spatial_niche_analysis.py](#run_spatial_niche_analysispy)
  - [run_de_analysis.py](#run_de_analysispy)
  - [colocation.py](#colocationpy)
  - [inspect_data_structure.py](#inspect_data_structurepy)
  - [tcga_survival_analysis.R](#tcga_survival_analysisr)
  - [requirements_r.txt](#requirements_rtxt)
- [各步骤运行命令](#各步骤运行命令)
- [当前代码中的不合理之处与改进建议](#当前代码中的不合理之处与改进建议)
- [输出文件汇总](#输出文件汇总)

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
scRNA-seq (GSE149614)          Visium 空间切片 (CHC20/CHC23)
        |                                  |
   [pre.py / preprocessing.py]             |
   数据预处理 & 质控                        质控 & 归一化 & 聚类
        |                                  |
        +---------> Cell2location <---------+
                  [run_preprocessing.py]
                  RegressionModel (scRNA) → 细胞类型参考签名
                  Cell2location (Visium) → 每个 spot 的细胞类型丰度
                         |
              [run_spatial_niche_analysis.py]
              构建免疫抑制空间生态位评分 & 基因签名
                         |
              [run_de_analysis.py] (可选 Wilcoxon 交叉验证)
                         |
              [tcga_survival_analysis.R]
              ssGSEA 投影 → KM + Cox 生存分析
```

---

## 文件结构说明

### `pre.py`

**功能定位**：早期数据预处理辅助脚本（前期探索性工具），主要用于从原始 txt 格式数据构建 scRNA-seq 参考数据，以及加载和合并多个 Visium 切片。

**主要函数**：

| 函数名 | 功能 |
|---|---|
| `build_scrna_reference()` | 读取 GSE149614 的 count 矩阵（txt）和 metadata，转为 AnnData 对象并保存为 h5ad |
| `_force_gene_id_index()` | 将 AnnData 的基因索引（var_names）强制切换为 Ensembl ID 风格，并去除版本号后缀 |
| `_load_visium_slice()` | 加载单个 Visium 切片，并将 gene_ids 列设为基因索引 |
| `load_and_merge_visium()` | 分别加载 CHC20 和 CHC23 两张 Visium 切片，拼接为单一 AnnData |
| `_align_shared_gene_ids()` | 求 scRNA 与 Visium 在 Ensembl ID 层面的交集，统一基因列表 |
| `diagnose_alignment()` | 打印详细的索引对齐诊断信息，方便调试两套数据的基因命名体系是否一致 |
| `_guess_id_style()` | 粗略判断基因索引是 SYMBOL 风格还是 Ensembl/gene_id 风格 |
| `_print_adata_snapshot()` | 打印 AnnData 的关键结构信息（形状、列名、前几行样本内容） |

**设计特点**：专门处理了 Ensembl ID 与 Gene Symbol 不一致的问题（空间组默认以 Symbol 为索引，scRNA-seq 原始数据以 Symbol 为主），统一采用 Ensembl ID 进行对齐，确保 Cell2location 建模的基因索引一致性。

---

### `preprocessing.py`

**功能定位**：核心预处理函数库，被 `run_preprocessing.py` 调用，封装了从 h5ad 加载数据到 Cell2location 建模所需的所有中间步骤。

**主要函数**：

| 函数名 | 功能 |
|---|---|
| `check_cell2location_available()` | 检查并导入 cell2location 包中的 RegressionModel 和 Cell2location 类 |
| `load_scrna_h5ad()` | 读取 scRNA h5ad 文件；过滤低质量细胞（总表达 < 200）和低表达基因（总表达 < 10）；统计各细胞类型数量 |
| `load_visium()` | 用 squidpy 读取 Visium 数据，计算 QC 指标（mt%），过滤 spot（总count 5000–35000，mt% < 20%，基因至少在 10 个 spot 表达），归一化、对数转换、高变基因选取、PCA、UMAP、Leiden 聚类 |
| `align_shared_genes()` | 对齐 scRNA 与 Visium 的共同基因（基于 Gene Symbol） |
| `subset_t_cells()` | 提取 T/NK 细胞子集，用于可视化确定 Treg 簇 |
| `assign_treg_label()` | 在新注释列中，将指定 T/NK 簇（默认为簇 9）重命名为 Treg，生成 `final_celltype` 列 |
| `setup_and_train_regression_model()` | 配置并训练 RegressionModel，从 scRNA 数据估算各细胞类型的参考表达签名（100 epochs） |
| `export_signatures()` | 导出 RegressionModel 的后验参考签名（写回 adata_sc.varm） |
| `extract_cell_state_df()` | 从 `adata_sc.varm["means_per_cluster_mu_fg"]` 中提取为标准 DataFrame |
| `sanitize_cell_state_df()` | 清理签名矩阵中的 inf/nan/负值，确保送入 Cell2location 的矩阵有效 |
| `compute_spot_cell_proportion()` | 将 Cell2location 输出的细胞丰度归一化为比例，生成 spot × 细胞类型的比例表 |

**设计特点**：所有功能高度模块化，每个函数职责单一，便于在主脚本中灵活组合调用。

---

### `run_preprocessing.py`

**功能定位**：整个分析的**第一步主运行脚本**，串联调用 `preprocessing.py` 中的功能模块，完成从原始数据到 Cell2location 反卷积结果的完整流程。

**主要流程**：

1. 读取 scRNA h5ad 参考数据（已预处理）
2. 读取 CHC20 Visium 空间切片（仅单张）
3. 对齐两套数据的共同基因
4. 可视化 T/NK 细胞中关键 Treg 标志物（CD3D、CD4、FOXP3、IL2RA），以点图辅助手动确认 Treg 簇
5. 将簇 9 的 T/NK 细胞重新注释为 Treg，生成 `final_celltype` 列
6. 训练 RegressionModel（100 epochs），获得细胞类型参考签名
7. 保存训练损失曲线图（`regression_training_history.png`）
8. 提取并清理签名矩阵（cell_state_df）
9. 以原始 count 矩阵运行 Cell2location 空间建模（1000 epochs，N_cells_per_location=30）
10. 保存空间建模训练曲线（`spatial_mapping_training_history.png`）
11. 导出后验细胞丰度，使用 `q05_cell_abundance_w_sf` 作为默认丰度表
12. 计算每个 spot 的细胞类型比例（`spot_cell_proportion.csv`）
13. 基于固定阈值（Hepatocyte > 0.3 且 Treg > 0.1）打简单共定位标签（`coloc` 列）
14. 保存 adata_sc_post.h5ad、adata_vis_post.h5ad、shared_genes.txt

**输入数据**：
- `data/scRNA_reference.h5ad`（由 `pre.py` 的 `build_scrna_reference()` 生成，或直接提供）
- `data/CHC20_Visium/`（10x Space Ranger 输出目录）

**主要输出**：
- `results/adata_sc_post.h5ad`
- `results/adata_vis_post.h5ad`
- `results/spot_cell_proportion.csv`
- `results/spot_with_coloc_label.csv`
- `results/regression_training_history.png`
- `results/spatial_mapping_training_history.png`

---

### `run_spatial_niche_analysis.py`

**功能定位**：整个分析的**第二步主运行脚本**，基于 Cell2location 反卷积得到的细胞丰度，定义免疫抑制空间生态位评分，并提取生态位特异性基因签名。

**免疫抑制相关基因面板（IMMUNOSUPPRESSIVE_GENES）**：
`FOXP3`, `IL2RA`, `CTLA4`, `TIGIT`, `LAG3`, `PDCD1`, `HAVCR2`, `TGFB1`, `IL10`, `CXCL12`, `CCL22`, `TNFRSF18`, `TNFRSF4`, `IKZF2`

**主要计算逻辑**：

| 指标 | 计算方式 |
|---|---|
| `hep_high` | Hepatocyte 比例 ≥ 75% 分位数的 spot，视为肿瘤实质区高丰度 spot |
| `distance_to_hep_high` | 每个 spot 到最近 hep_high spot 的欧氏距离（用 cKDTree 计算） |
| `Treg_like_score` | Z-score(Treg) + Z-score(免疫抑制基因模块得分) |
| `immune_stroma_score` | Z-score(Treg) + Z-score(T/NK) + Z-score(Myeloid) + Z-score(Fibroblast) |
| `immunosuppressive_niche_score` | Z-score(Treg) + Z-score(Myeloid) + Z-score(Fibroblast) + Z-score(免疫抑制基因模块) + Z-score(邻域 Hepatocyte 均值) |
| `niche_high` | immunosuppressive_niche_score ≥ 80% 分位数的 spot |
| `spatial_region` | 基于 hep_high 与 immune_stroma_score 划分为：tumor_core / tumor_edge / stroma_immune / other |

**空间邻域计算**：以每个 spot 到其最近邻居的中位距离 × 1.25 为半径，用 KD-tree 检索邻域 spot，计算邻域均值特征。

**基因签名提取**：对 niche_high 与 niche_low 两组 spot 分别计算 log2 fold change，按 log2FC 排序后输出前 80 个基因。

**输出文件**：
- `results/spatial_niche/spatial_niche_scores.csv`（每个 spot 的全部评分指标）
- `results/spatial_niche/spatial_niche_parameters.csv`（所用阈值和参数记录）
- `results/spatial_niche/immunosuppressive_niche_signature_genes.txt`
- `results/spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv`
- `results/spatial_signature_genes.txt`（供 TCGA 生存分析使用）
- `results/spatial_niche/plots/` 目录下 10 张可视化图

---

### `run_de_analysis.py`

**功能定位**：**可选的第三步**，对 Cell2location 反卷积后的 Visium 数据，使用 Scanpy 的 Wilcoxon 秩和检验进行差异表达分析，提供独立于 log2FC 排序的签名基因交叉验证。

**主要流程**：

1. 读取 `adata_vis_post.h5ad`（Cell2location 结果）
2. 读取共定位标签文件（CSV），对齐到 adata 的 obs_names
3. 若指定 layer（如 `log1p`）存在则直接使用；否则自动进行 normalize_total + log1p 标准化
4. 运行 `sc.tl.rank_genes_groups()` 进行 Wilcoxon DE 检验
5. 提取共定位阳性组（coloc=1 或 niche_high=True）的前 N 个差异表达基因
6. 输出到签名基因文件（`spatial_signature_genes.txt`）

**灵活性**：通过命令行参数支持不同的共定位标签文件（`--coloc`）和标签列（`--coloc-column`），可与 `run_preprocessing.py` 产生的简单阈值标签或 `run_spatial_niche_analysis.py` 产生的生态位高评分标签配合使用。

---

### `colocation.py`

**功能定位**：一个**独立的共定位辅助分析脚本**，提供了基于逻辑回归的软性共定位评分方法，作为 `run_preprocessing.py` 中硬阈值方法的替代方案。

**核心函数 `add_coloc_logistic()`**：

1. 以"肿瘤细胞比例 + Treg 比例之和"超过某一阈值（支持自动分位数估计）为临时正样本目标
2. 用逻辑回归模型学习两种细胞比例与共定位可能性的关系
3. 输出每个 spot 的共定位概率分数（`coloc_score`）
4. 基于概率分数的分位数（默认 80%）打硬标签（`coloc`）

**与主流程的关系**：该脚本独立于主流程之外，目前未被 `run_preprocessing.py` 或 `run_spatial_niche_analysis.py` 直接调用，可作为可选的替代共定位定义方法单独运行。

---

### `inspect_data_structure.py`

**功能定位**：**调试辅助工具**，用于在分析早期快速检查 AnnData 对象的数据结构，帮助开发者了解数据格式，确认列名、obs/var 内容等信息。

**主要函数**：

| 函数名 | 功能 |
|---|---|
| `print_structure()` | 打印 AnnData 的 shape、obs 列、var 列、uns 键、layers 键 |
| `print_query()` | 按需打印 obs/var 的前 N 行内容，支持 sc_obs、sc_var、sp_obs、sp_var 四种查询模式 |

**使用方式**（命令行）：

```bash
python inspect_data_structure.py \
    --sc-h5ad results/adata_sc_post.h5ad \
    --spatial-h5ad results/adata_vis_post.h5ad \
    --query sc_obs --head 10
```

---

### `tcga_survival_analysis.R`

**功能定位**：**第四步和第五步**，将空间生态位基因签名投影到 TCGA-LIHC（肝细胞癌）队列的 bulk RNA-seq 数据，进行单样本 GSEA（ssGSEA）评分和 Kaplan-Meier + Cox 生存分析。

**主要流程**：

1. **参数解析**：支持命令行参数传入数据集、签名文件路径、输出目录等（兼容 Windows / Linux / macOS 路径）
2. **GDC 数据下载**：通过 `TCGAbiolinks` 从 GDC 下载 TCGA-LIHC 的 RNA-seq 数据（HTSeq-FPKM），带重试逻辑（默认重试 3 次），支持 SSL 证书验证控制，支持缓存 RDS 避免重复下载
3. **表达数据预处理**：处理重复基因名（按均值聚合），log2(FPKM+1) 变换
4. **ssGSEA 评分**（Step 4）：使用 `GSVA::gsva()` 以 ssGSEA 方法计算每个样本的空间签名评分，按患者 ID 聚合（取均值）
5. **生存分析**（Step 5）：
   - 下载临床数据，提取 `time`（days_to_death 或 days_to_last_follow_up）和 `event`（死亡事件）
   - 自动探测临床表中的年龄（age_at_diagnosis）、肿瘤分期（ajcc_pathologic_stage）等协变量
   - 拟合多变量 Cox 比例风险模型，输出结果到文本文件
   - 按评分中位数分为 High/Low 组，绘制 Kaplan-Meier 生存曲线（含 p 值和风险表）
   - 绘制 Cox 森林图（Forest Plot）

**输出文件**：
- `results/tcga_signature_score.csv`（每个患者的 ssGSEA 评分）
- `results/tcga_signature_survival.csv`（评分与临床数据合并表）
- `results/cox_results.txt`（Cox 模型摘要）
- `results/km_plot.png`（KM 生存曲线图）
- `results/cox_forest_plot.png`（Cox 森林图）

---

### `requirements_r.txt`

R 包依赖列表：

```
TCGAbiolinks
GSVA
survival
survminer
SummarizedExperiment
```

---

## 各步骤运行命令

> 以下命令以 Linux 服务器环境为例，Python 环境路径和数据路径请按实际修改。

### Step 1：Cell2location 细胞类型反卷积（含 scRNA 预处理）

```bash
/mnt/hdd/private/gyw/.virtualenvs/Single_cell_train/bin/python code/run_preprocessing.py
```

若需要先从 txt 原始数据构建 scRNA 参考数据，运行：

```bash
python code/pre.py
```

### Step 2：空间免疫抑制生态位分析

```bash
/mnt/hdd/private/gyw/.virtualenvs/Single_cell_train/bin/python code/run_spatial_niche_analysis.py
```

支持的可选参数示例：

```bash
python code/run_spatial_niche_analysis.py \
    --adata results/adata_vis_post.h5ad \
    --hep-high-quantile 0.75 \
    --niche-high-quantile 0.80 \
    --top-niche-genes 80
```

### Step 3（可选）：Wilcoxon 签名基因交叉验证

```bash
/mnt/hdd/private/gyw/.virtualenvs/Single_cell_train/bin/python code/run_de_analysis.py \
    --coloc results/spatial_niche/spatial_niche_scores.csv \
    --coloc-column niche_high \
    --out results/spatial_niche/wilcoxon_niche_signature_genes.txt \
    --top-n 80
```

### Step 3（可选替代）：基于逻辑回归的共定位评分

```bash
python code/colocation.py \
    --input results/spot_cell_proportion.csv \
    --output results/spot_with_coloc_label.csv \
    --tumor-col Hepatocyte --treg-col Treg \
    --label-quantile 0.8
```

### Step 4 & 5：TCGA 生存分析（R）

```bash
Rscript code/tcga_survival_analysis.R \
    --signature results/spatial_signature_genes.txt \
    --out-dir results
```

---

## 当前代码中的不合理之处与改进建议

### 1. 硬编码路径问题（高优先级）

**问题**：`pre.py`、`run_preprocessing.py`、`colocation.py` 中存在大量硬编码的绝对路径，且不同文件中指向不同机器的路径：
- `pre.py` 中 `ROOT_DIR = Path("/home/gyw/R/Python/Single_cell_train")`
- `run_preprocessing.py` 中 `ROOT_DIR = Path("/mnt/hdd/private/gyw/Code/txbb/Single_cell_train")`
- `colocation.py` 中 `--input default="/home/gyw/R/Python/Single_cell_train/results/spot_cell_proportion.csv"`
- `tcga_survival_analysis.R` 中 `default_root = "E:\\Research!!\\Codes\\Bioinformation_train\\Single-cell&spatial_transcriptomics"`

这些路径与个人开发机强绑定，导致代码无法直接在其他环境运行。

**建议**：统一通过脚本自身位置（`Path(__file__).parent`）动态推断项目根目录，或引入 `.env` 配置文件管理路径，参照 `run_de_analysis.py` 和 `run_spatial_niche_analysis.py` 的良好实践，所有脚本应统一使用 `Path(__file__).resolve().parents[1]` 的方式计算路径。

---

### 2. 两套预处理模块并存（高优先级）

**问题**：`pre.py` 和 `preprocessing.py` 存在大量功能重叠，且部分功能存在不一致：
- `pre.py` 使用 Ensembl ID 对齐（更严格），`preprocessing.py` 使用 Gene Symbol 对齐（`align_shared_genes()` 直接对 `var_names` 取交集，而 `run_preprocessing.py` 读取的是已经是 Symbol 索引的 h5ad）
- `pre.py` 中 `_load_visium_slice()` 使用 `sc.read_visium()`，`preprocessing.py` 中 `load_visium()` 使用 `sq.read.visium()`，两者接口不同
- `pre.py` 有专门的 `diagnose_alignment()` 诊断工具，但 `preprocessing.py` 中没有对应检查

**建议**：合并为单一预处理模块，明确以 Ensembl ID 为统一基因索引标准，删除冗余代码。

---

### 3. RegressionModel 训练轮次过少（中优先级）

**问题**：`run_preprocessing.py` 中 RegressionModel 训练 `max_epochs=100`，而 `pre.py` 中使用 `max_epochs=250`，官方推荐通常为 200–300 epochs，且没有 Early Stopping。100 epochs 可能导致模型未充分收敛，影响细胞类型签名的准确性。

**建议**：将 `max_epochs` 统一调整为 250，或添加 early stopping 回调，并在训练后检查损失曲线是否趋于平稳。

---

### 4. Cell2location 仅处理单张切片（中优先级）

**问题**：`run_preprocessing.py` 中 `load_visium(PATH_CHC20, sample_name="CHC20")` 只加载了 CHC20 一张切片，而 `pre.py` 中已提供了完整的 `load_and_merge_visium()` 函数可加载 CHC20+CHC23 两张切片并合并。单样本分析限制了分析的普适性和统计功效。

**建议**：在主流程中使用 `pre.py` 的多切片合并逻辑，对两张 Visium 切片同时进行分析，并在后续空间分析中加入样本来源（`sample` 列）作为批次校正或分层分析的依据。

---

### 5. 简单硬阈值共定位标签的生物学局限性（中优先级）

**问题**：`run_preprocessing.py` 中的共定位标签定义为 `(Hepatocyte > 0.3) & (Treg > 0.1)` 的固定阈值布尔判断，在 README 中也已注明这是旧方法并不再推荐。然而这段代码仍然被运行并输出 `spot_with_coloc_label.csv`，但后续 `run_de_analysis.py` 默认读取该文件时，实际上在新流程中已被 `run_spatial_niche_analysis.py` 的 `niche_high` 列替代。两套标签系统并存容易造成混淆。

**建议**：移除 `run_preprocessing.py` 中 Step 2 的硬阈值共定位代码，或将其标记为仅用于快速预览的参考。统一推荐使用 `run_spatial_niche_analysis.py` 产生的 `niche_high` 标签作为权威共定位定义。

---

### 6. `q05` 与 `means` 丰度选择缺乏明确说明（中优先级）

**问题**：`run_preprocessing.py` 第 102 行有注释 `# TODO: 这里太保守`，将 `q05_cell_abundance_w_sf`（后验 5% 分位数，非常保守）赋值给 `cell_abundance`，而 `run_spatial_niche_analysis.py` 的参数默认值为 `means_cell_abundance_w_sf`（后验均值）。两者在下游比例计算中会产生不同结果，且存在不一致性。

**建议**：统一使用 `means_cell_abundance_w_sf` 作为默认丰度指标（后验均值），q05 可保留作为保守估计的备用选项，并在参数文档中明确说明两者的适用场景。

---

### 7. `_plot_hep_treg()` 中 Treg 分位数参考线使用 75% 而非参数中的 niche 分位数（低优先级）

**问题**：`run_spatial_niche_analysis.py` 中 `_plot_hep_treg()` 函数内，参考线使用 `df["Treg"].quantile(0.75)` 硬编码为 75% 分位数，而 `hep_high_quantile` 参数的默认值也是 0.75，但 `niche_high_quantile` 参数为 0.80。图中的分位数参考线没有与实际用于定义 `hep_high` 的分位数阈值（`hep_thr`）一致。

**建议**：将 `_plot_hep_treg()` 中的参考线改为传入实际使用的 `hep_thr` 和 `treg_thr`（而非重新计算分位数），确保图与分析参数一致。

---

### 8. `colocation.py` 的逻辑回归方法存在数据泄露风险（低优先级）

**问题**：`add_coloc_logistic()` 函数以两种细胞比例之和构建临时标签 `y`，再用同样的两个特征训练逻辑回归模型预测 `y`。训练目标直接由特征线性组合而来，模型学到的本质上是同一线性边界，与直接设阈值等价，不具有真正的机器学习泛化价值。

**建议**：若要使用模型方法定义共定位，应结合更多生物学特征（基因表达模块评分、空间邻域特征等），或改用无监督方法（如 GMM、DBSCAN）直接对丰度空间进行聚类。

---

### 9. `run_de_analysis.py` 在数据未归一化时进行 DE 可能出错（低优先级）

**问题**：当 `adata_vis_post.h5ad` 中不存在 `log1p` layer 时，脚本会对 `adata.X` 直接进行 `normalize_total + log1p`，但此时 `adata.X` 可能已经是 Cell2location 处理后的数据（可能包含非整数值或已归一化数据），再次归一化会引入错误。

**建议**：明确要求调用者传入一个存储有原始 count 或已正确归一化的 layer，或在脚本开头检查 `adata.X` 的数据类型和范围，给出明确的警告信息。

---

### 10. 缺少统一的 Python 环境依赖文件（低优先级）

**问题**：项目提供了 R 包的 `requirements_r.txt`，但缺少 Python 包的 `requirements.txt` 或 `environment.yml`，无法一键还原 Python 环境。

**建议**：添加 `requirements.txt`（或 `environment.yml`），包含：
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
```

---

## 输出文件汇总

| 步骤 | 输出文件 | 说明 |
|---|---|---|
| Step 1 | `results/adata_sc_post.h5ad` | 训练后的 scRNA AnnData |
| Step 1 | `results/adata_vis_post.h5ad` | Cell2location 反卷积后的 Visium AnnData |
| Step 1 | `results/spot_cell_proportion.csv` | 每个 spot 的细胞类型比例表 |
| Step 1 | `results/regression_training_history.png` | RegressionModel 训练损失曲线 |
| Step 1 | `results/spatial_mapping_training_history.png` | Cell2location 训练损失曲线 |
| Step 1 | `results/t_cell_dotplot.png` | T/NK 簇标志基因点图（确认 Treg 簇） |
| Step 1 | `results/shared_genes.txt` | scRNA 与 Visium 共同基因列表 |
| Step 2 | `results/spatial_niche/spatial_niche_scores.csv` | 每个 spot 的全部生态位评分 |
| Step 2 | `results/spatial_niche/spatial_niche_parameters.csv` | 分析所用参数与阈值记录 |
| Step 2 | `results/spatial_niche/immunosuppressive_niche_signature_genes.txt` | 免疫抑制生态位签名基因（用于 TCGA） |
| Step 2 | `results/spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv` | 带 log2FC 的完整排序签名表 |
| Step 2 | `results/spatial_signature_genes.txt` | TCGA 分析直接使用的签名基因文件 |
| Step 2 | `results/spatial_niche/plots/*.png` | 10 张空间可视化图 |
| Step 3 | `results/spatial_niche/wilcoxon_niche_signature_genes.txt` | Wilcoxon DE 验证签名基因（可选） |
| Step 4-5 | `results/tcga_signature_score.csv` | 每个 TCGA 患者的 ssGSEA 评分 |
| Step 4-5 | `results/tcga_signature_survival.csv` | 评分与生存数据合并表 |
| Step 4-5 | `results/cox_results.txt` | Cox 多变量回归结果 |
| Step 4-5 | `results/km_plot.png` | Kaplan-Meier 生存曲线图 |
| Step 4-5 | `results/cox_forest_plot.png` | Cox 森林图 |
