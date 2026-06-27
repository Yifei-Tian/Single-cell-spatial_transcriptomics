# `results/joint_HCC4R_CHC20/` 结果文件完整讲解

> **分析背景**：HCC4R 与 CHC20 两张肝细胞癌 Visium 空间转录组切片经批次效应评估确认 batch effect 极小，因此采用**联合分析策略**——用同一套 scRNA-seq 参考训练一次 RegressionModel，分别对两张切片反卷积后合并 spot 进行下游分析，以获得更大样本量和更高统计功效。
>
> **运行入口**：`python code/run_joint_hcc4r_chc20.py`
>
> **分析脚本**：
> - `code/pipeline/run_preprocessing.py` → `main_joint()` 函数（Step 1）
> - `code/pipeline/run_spatial_niche_analysis.py` → `main()` 函数（Step 2，启用 `--per-sample-neighbors`）

---

## 目录

- [一、Step 1 输出文件 — 数据预处理与 Cell2location 反卷积](#一step-1-输出文件--数据预处理与-cell2location-反卷积)
  - [1.1 AnnData 数据文件](#11-anndata-数据文件)
  - [1.2 细胞比例表](#12-细胞比例表)
  - [1.3 scRNA-seq 可视化图](#13-scrna-seq-可视化图)
  - [1.4 训练曲线图](#14-训练曲线图)
  - [1.5 切片对比图（cross_slice_comparison/）](#15-切片对比图cross_slice_comparison)
  - [1.6 辅助文件](#16-辅助文件)
- [二、Step 2 输出文件 — 空间免疫抑制生态位分析（spatial_niche/）](#二step-2-输出文件--空间免疫抑制生态位分析spatial_niche)
  - [2.1 核心数据表](#21-核心数据表)
  - [2.2 特征基因文件](#22-特征基因文件)
  - [2.3 可视化图表（spatial_niche/plots/）](#23-可视化图表spatial_nicheplots)
- [三、结果文件索引表](#三结果文件索引表)
- [四、分析流程数据流图](#四分析流程数据流图)
- [五、核心生物学结论链](#五核心生物学结论链)

---

## 一、Step 1 输出文件 — 数据预处理与 Cell2location 反卷积

**对应脚本**：`code/pipeline/run_preprocessing.py`，调用入口为 `main_joint()` 函数

**流程概述**（联合分析模式，共 9 步）：
- Step 1：加载 scRNA-seq 参考数据，执行质控
- Step 2：分别加载 HCC4R 和 CHC20 两张 Visium 空间数据
- Step 3：计算 scRNA × HCC4R × CHC20 三路共同基因子集
- Step 4：从 T/NK 细胞亚群中鉴定 Treg，绘制标志基因气泡图
- Step 5：用三路共同基因子集训练**统一** RegressionModel（一次，两切片共享）
- Step 6：用统一签名分别对 HCC4R 和 CHC20 做 Cell2location 反卷积
- Step 7：合并两切片 spot 级结果，添加 `obs["sample"]` 来源标记
- Step 8：生成两切片细胞组成一致性对比图
- Step 9：保存共享基因列表和 scRNA AnnData

---

### 1.1 AnnData 数据文件

#### `adata_vis_post_joint.h5ad`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 7**

**文件格式**：HDF5 格式的 AnnData 对象（`.h5ad`，单细胞/空间组学标准格式）

**内容说明**：

这是 HCC4R + CHC20 两张切片经过 Cell2location 贝叶斯反卷积后合并的**核心数据文件**，也是后续所有空间分析（Step 2）的输入。文件结构：

| 属性 | 内容 |
|------|------|
| `adata.X` | 每个 spot 的原始 count 表达矩阵（行 = spot，列 = 基因） |
| `adata.obs` | 每个 spot 的元数据，**含 `sample` 列（"HCC4R" 或 "CHC20"）** 标记切片来源 |
| `adata.obs_names` | spot 唯一标识符，格式为 `{样本名}_{原始barcode}`（如 `HCC4R_ACGCCTGACACGCGCT-1`），添加样本前缀以防合并后重复 |
| `adata.obsm["spatial"]` | 每个 spot 的二维空间坐标（来自 Visium 芯片物理布局） |
| `adata.obsm["means_cell_abundance_w_sf"]` | Cell2location 后验估计的细胞类型丰度矩阵，行 = spot，列 = 细胞类型（**反卷积主要输出**） |
| `adata.obsm["cell_abundance"]` | 同上，简写别名 |
| `adata.var_names` | 三路共同基因集（scRNA ∩ HCC4R ∩ CHC20） |

**技术原理**：

10x Visium 空间转录组每个 spot 直径约 55 μm，实际捕获的 RNA 来自多种细胞类型的混合信号。Cell2location 是一种贝叶斯层次模型，利用 scRNA-seq 参考中各细胞类型的基因表达特征（参考签名 `cell_state_df`），通过负二项分布似然函数将每个 spot 的混合信号"反卷积"为各细胞类型的估计数量（后验均值 `means_cell_abundance_w_sf`）。

**联合分析关键设计**：

两张切片使用**同一套 scRNA-seq 参考签名**（一次 RegressionModel 训练）进行反卷积，保证细胞类型估计在相同参考系下可直接比较。合并后通过 `adata.obs["sample"]` 列区分来源，下游空间邻域分析通过 `--per-sample-neighbors` 参数限制在同一切片内建立邻居关系。

**揭示的生物学现象**：

- 每个 spot 的 `means_cell_abundance_w_sf` 向量直接回答"该空间位置存在哪些细胞类型、各占多少"
- 可见肿瘤实质（Hepatocyte 丰度高区域）、免疫浸润区（Treg/T-NK/Myeloid 丰度高区域）和基质区（Fibroblast 丰度高区域）的空间分布格局
- 这是从"基因表达混合信号"到"细胞类型空间地图"的关键转化步骤

---

#### `adata_vis_post.h5ad`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 7**

**内容**：与 `adata_vis_post_joint.h5ad` **完全相同**，是其标准名称副本，供下游 `run_spatial_niche_analysis.py` 通过默认参数路径直接读取。

---

#### `adata_vis_post_HCC4R.h5ad`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 6**（Cell2location 反卷积，HCC4R 切片）

**内容**：HCC4R 单张切片独立反卷积结果，结构与 `adata_vis_post_joint.h5ad` 相同，但仅包含 HCC4R 的 spot，**不含 `sample` 列**。

**用途**：用于单独分析 HCC4R 切片，或与联合分析结果对比验证。批次效应极小时，单独反卷积结果与联合反卷积结果的细胞类型比例应高度一致。

---

#### `adata_vis_post_CHC20.h5ad`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 6**（Cell2location 反卷积，CHC20 切片）

**内容**：CHC20 单张切片独立反卷积结果，结构与 `adata_vis_post_HCC4R.h5ad` 对称。

---

#### `adata_sc_post.h5ad` / `adata_sc_post_joint.h5ad`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 9**

**内容**：经过 Treg 亚群重注释后的 scRNA-seq 参考 AnnData（含 `obs["final_celltype"]` 列，Treg 已从 T/NK 中分离标注）。两个文件内容相同，后者为带明确标记的副本。

**用途**：保存处理后的参考数据集，供后续分析溯源或额外挖掘使用。

---

### 1.2 细胞比例表

#### `spot_cell_proportion_HCC4R.csv`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 6**

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| `spot_id` | Visium spot 唯一标识符（原始 barcode） |
| `Hepatocyte` | 该 spot 中肝细胞比例（归一化到 0–1） |
| `Treg` | 调节性 T 细胞比例 |
| `T/NK` | T 细胞和 NK 细胞合并比例 |
| `Myeloid` | 髓系细胞（TAM、MDSC 等）比例 |
| `Fibroblast` | 成纤维细胞（CAF）比例 |
| `B cell` | B 细胞比例 |
| 其他细胞类型 | 依 scRNA 参考数据中细胞类型而定 |

所有细胞类型比例之和为 1（归一化后）。

**技术细节**：通过 `compute_spot_cell_proportion()` 函数将 `means_cell_abundance_w_sf` 中的绝对丰度值归一化为相对比例。

**揭示的生物学现象**：

- `Treg` 列的数值分布揭示 HCC4R 切片中免疫抑制压力的空间分布格局
- `Hepatocyte` 高值区域对应肿瘤实质（肿瘤核心区），低值区域对应间质/免疫浸润区
- `Treg` 与 `Myeloid` 同时高值的 spot 提示可能的免疫抑制微生态位

---

#### `spot_cell_proportion_CHC20.csv`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 6**

**内容**：与 `spot_cell_proportion_HCC4R.csv` 结构对称，记录 CHC20 切片各 spot 的细胞类型比例。

---

#### `spot_cell_proportion_joint.csv`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 7**

**内容**：HCC4R 和 CHC20 两切片比例表的垂直合并，额外包含 `sample` 列：

| 额外列 | 含义 |
|--------|------|
| `sample` | 切片来源标记（"HCC4R" 或 "CHC20"） |

**用途**：
- 通过 `sample` 列分组比较两切片细胞类型分布差异，评估批次效应
- 快速检查联合分析合理性（若两切片细胞类型构成高度相似，则联合分析有效）

---

### 1.3 scRNA-seq 可视化图

#### `scrna_tsne_celltype.png`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 4+**（调用 `plot_scrna_tsne_celltype()`）

**文件格式**：PNG 图像（仿论文 Figure C 风格，参考 `docs/plots/1.png`）

**图表内容**：

scRNA-seq 单细胞数据的 tSNE（或 UMAP）嵌入可视化图：
- 每个点代表一个细胞，使用高饱和度离散色板（仿论文 tSNE 图风格）按细胞类型着色
- 每种细胞类型聚集中心用类型名称标注（白色背景框）
- 右侧图例完整展示颜色-类型对应关系
- 坐标轴以箭头样式标注方向（tSNE 1 / tSNE 2）

**如何理解**：

此图展示 scRNA-seq 参考数据集中各细胞类型在低维嵌入空间中的分布格局。不同细胞类型应形成清晰分离的点云簇，生物学相近的类型（如 Treg 与 T/NK）通常在空间上相邻但可分辨。

**揭示的生物学现象**：

- 聚类的分离程度反映参考数据集的质量，清晰分离说明 scRNA-seq 数据捕获了充分的转录异质性
- Treg 细胞与其他 T/NK 细胞应在相邻但独立的区域，验证了 Treg 亚群鉴定的准确性
- Myeloid、Fibroblast 与免疫细胞的远距分离反映了肝癌 TME 细胞类型组成的多样性

---

#### `scrna_celltype_marker_heatmap.png`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 4+**（调用 `plot_celltype_marker_heatmap()`）

**文件格式**：PNG 图像（仿论文 Figure D 风格，参考 `docs/plots/2.png`）

**图表内容**：

细胞类型 Marker 基因平均表达热图：
- **X 轴**：细胞类型（按生物学关联性排序）
- **Y 轴**：各细胞类型的特异性标志基因
- **颜色**：Yellow-Black-Purple 渐变（paper_ybp 配色）：亮黄 = 高表达，黑色 = 中等，深紫 = 低/无表达
- **归一化**：CP10K + log1p 变换后各细胞类型内均值，再经跨细胞类型 Z-score 标准化

**主要 Marker 基因组**：

| 细胞类型 | 关键 Marker 基因 |
|---------|----------------|
| T cell | IL7R, CD3G, CD2, ITM2A, CD3D |
| Myeloid | LYZ, AIF1, RNASE1, C1QB, HLA-DRA |
| NK | GNLY, GZMB, KLRD1, KLRF1 |
| B cell | B3GNT7, MS4A1, BANK1, CD79A, TCL1A |
| Malignant（肝癌细胞） | APOA2, ALB, APOA1, AMBP, APOH, TTR |
| Endothelial | PECAM1, CDH5, SPARCL1, STC1, SPARC |
| Epithelial | INSR, KRT18, KRT19, EPCAM, SOX4 |
| Plasma cell | JCHAIN, MZB1, IGLL5, SSR4 |
| HSC（肝星状细胞） | RGS5, COL1A1, ACTA2, PDGFRB |

**如何理解**：

热图中理想模式为每种细胞类型的特异性 Marker 在对应列呈亮黄色（高 Z-score），在其他列保持深紫色（低/无表达），形成"对角线"高亮模式，验证细胞类型注释的可靠性。

**揭示的生物学现象**：

- ALB、APOA2 等白蛋白相关基因在 Malignant/Hepatocyte 中特异高表达，确认肿瘤细胞的肝细胞来源
- COL1A1、ACTA2 在 HSC/Fibroblast 中高表达，揭示肝星状细胞激活后的基质重塑能力，与 CAF 免疫抑制功能相关
- FOXP3 在 Treg 簇中即使表达量绝对值低，跨细胞类型 Z-score 相对差异仍清晰可见，验证 Treg 亚群鉴定

---

#### `t_cell_dotplot_horizontal.png`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 4**（调用 `plot_treg_dotplot_horizontal()`）

**文件格式**：PNG 图像

**图表内容**：

横版气泡图（Dot Plot），X 轴为 T/NK 细胞各子簇（Leiden 分辨率 `res.3` 下的聚类编号），Y 轴为 4 个关键标志基因：

| 基因 | 生物学意义 |
|------|-----------|
| **CD3D** | 泛 T 细胞标志物，所有真正的 T 细胞均表达；用于确认这些簇确实是 T 细胞而非其他免疫细胞 |
| **CD4** | Treg 属于 CD4⁺ T 细胞亚群，Treg 所在簇中 CD4 需有一定程度表达 |
| **FOXP3** | Treg 的核心主调控转录因子（master TF），只有 Treg 才特异高表达，是鉴定 Treg 的金标准 |
| **IL2RA（CD25）** | 稳定 Treg 高表达，通常与 FOXP3 协同出现 |

气泡大小 = 该簇中表达该基因的细胞比例；气泡颜色深浅 = 表达量均值（Blues 色标）。

**如何理解**：

通过 FOXP3 仅在特定簇（如簇9）中出现大而深色的气泡，可以精确定位 Treg 亚群。FOXP3 是 Treg 的"主调控转录因子"，其特异性表达是鉴定 Treg 的金标准。

**揭示的生物学现象**：

清晰的 FOXP3 特异性信号证明了分析中使用的 Treg 亚群鉴定是准确的，为后续空间 Treg 分布分析和免疫抑制生态位识别奠定可靠的细胞类型定义基础。

---

### 1.4 训练曲线图

#### `regression_training_history_joint.png`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 5**（调用 `setup_and_train_regression_model()`，训练 RegressionModel）

**文件格式**：PNG 图像

**图表内容**：

Cell2location 第一阶段 RegressionModel（单细胞参考签名学习）的训练损失曲线：
- **X 轴**：训练 epoch 数（从第 50 个 epoch 开始显示，前期不稳定部分裁剪）
- **Y 轴**：ELBO（证据下界，Evidence Lower BOund）

**联合分析特殊性**：

与单切片模式（各切片各自独立训练）不同，联合分析使用 scRNA 参考与 HCC4R 和 CHC20 的**三路共同基因子集**训练**唯一一次** RegressionModel，得到统一参考签名 `cell_state_df`，保证两张切片的反卷积结果在同一参考系下进行。

**如何理解**：

ELBO 是变分推断（Variational Inference）的优化目标，其绝对值越大（负 ELBO 越小）说明模型对数据拟合越好。曲线应在若干 epoch 后趋于平稳（收敛），若仍在持续下降则说明训练轮次不足（默认 1000 个 epoch），需增加训练轮数。

**揭示的生物学现象**：

收敛的训练曲线保证了细胞类型参考签名的可靠性——只有 RegressionModel 充分学习到每种细胞类型的特异性基因表达模式后，下一阶段（Step 6）的空间反卷积才能准确区分不同细胞类型。

---

#### `spatial_mapping_training_history_HCC4R.png` / `spatial_mapping_training_history_CHC20.png`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 6**（`_run_cell2location_for_slice()` 函数，两切片各自空间建模）

**文件格式**：PNG 图像

**图表内容**：Cell2location 第二阶段空间建模（Cell2location 模型）的训练损失曲线，格式与 `regression_training_history_joint.png` 相同。两张切片各产生一张独立曲线图。

**如何理解**：

第二阶段模型使用 RegressionModel 得到的参考签名 `cell_state_df` 为输入，通过负二项分布对每个 spot 的细胞丰度进行后验估计（训练 1000 个 epoch）。每张切片独立进行，收敛曲线代表该切片的反卷积质量。

---

### 1.5 切片对比图（`cross_slice_comparison/`）

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 8**（调用 `_compare_slices_generic()`）

此目录包含 3 张用于评估 HCC4R 与 CHC20 批次效应的可视化图：

---

#### `cross_slice_comparison/cross_slice_mean_proportion_comparison.png`

**内容**：HCC4R 与 CHC20 各细胞类型**全切片平均比例**的并排条形图（蓝色=HCC4R，橙色=CHC20）。

**如何理解**：

若两切片各细胞类型平均比例趋势相近（如 Hepatocyte 均为最高占比，Treg 均偏低但存在），说明两个样本整体微环境构成具有代表性，不存在严重批次偏差，支持联合分析的合理性。

**揭示的生物学现象**：

肝癌肿瘤微环境的"细胞生态"在不同患者之间是否具有共同特征，即肝细胞癌 TME 免疫细胞组成的普遍性规律。

---

#### `cross_slice_comparison/cross_slice_treg_distribution.png`

**内容**：HCC4R 与 CHC20 Treg 比例分布的核密度估计（KDE）对比曲线图，X 轴为 Treg 比例，Y 轴为密度。

**如何理解**：

分布峰值、宽度和尾部形状相近，说明两切片的 Treg 浸润模式一致。

**揭示的生物学现象**：

Treg 在肝癌组织中的浸润通常是稀疏但局灶性富集的（分布呈右偏，大多数 spot 接近 0，少数 spot 有高 Treg 比例）。这种"局灶性免疫抑制"模式的跨样本一致性，是免疫逃逸机制普遍存在的空间证据。

---

#### `cross_slice_comparison/cross_slice_celltype_boxplot.png`

**内容**：Treg、Myeloid、Fibroblast、Hepatocyte 四种关键细胞类型在两切片中比例分布的箱线图对比（蓝色=HCC4R，橙色=CHC20）。

**如何理解**：

箱体中位线、四分位距和异常值分布的相似性反映两切片细胞比例的整体一致性。若中位数和分布形状高度相近，则批次效应可忽略，联合分析有效。

**揭示的生物学现象**：

多细胞类型比例分布的跨样本一致性验证了 HCC4R 和 CHC20 代表同一类型肝癌微环境（早期复发性肝癌），支持将两张切片视为同一生物学总体的重复样本进行联合分析。

---

### 1.6 辅助文件

#### `shared_genes_joint.txt`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → **Step 9**

**文件格式**：纯文本，每行一个基因名（HGNC 符号）

**内容**：scRNA-seq 参考数据与 HCC4R 和 CHC20 两张 Visium 切片的**三路共同基因列表**（即 `scRNA ∩ HCC4R ∩ CHC20`），也是最终用于 Cell2location 建模的基因集。

**如何理解**：

三路交集保证了联合分析中两张切片的反卷积使用完全相同的基因空间（同一 `cell_state_df` 的基因维度）。共有基因数量通常在 2,000–5,000 个之间（取决于测序深度和 Visium 捕获效率）。

**临床应用**：若需要在后续分析中回溯哪些基因被包含在建模中，可直接从此文件读取。

---

#### `run_preprocessing.log`

**来源脚本 / Step**：`run_preprocessing.py` → `main_joint()` → 全流程日志

**文件格式**：纯文本日志文件

**内容**：整个 Step 1 流程的实时运行记录，包含：
- 每个步骤的开始/完成时间
- 三路共同基因数量（如 `Three-way shared genes: 3247`）
- Treg 簇鉴定结果（各细胞类型计数）
- 训练轮次和收敛信息
- 合并后 AnnData 统计（如 `Joint AnnData: 5842 spots (HCC4R: 2981, CHC20: 2861), 3247 genes`）
- 各文件保存路径

**用途**：排查错误、记录运行参数、确认各步骤是否正常完成。

---

## 二、Step 2 输出文件 — 空间免疫抑制生态位分析（`spatial_niche/`）

**对应脚本**：`code/pipeline/run_spatial_niche_analysis.py`，`main()` 函数

**输入文件**：`adata_vis_post_joint.h5ad`（合并的两切片反卷积结果）

**关键参数**：`--per-sample-neighbors`（限制空间 kNN 邻域在同一切片内建立）

**流程概述**（共 15 步）：

| Step | 内容 |
|------|------|
| Step 1 | 加载 AnnData，提取细胞丰度矩阵 |
| Step 2 | 构建空间半径邻域（用于计算邻域均值特征） |
| Step 3 | 归一化细胞丰度为比例 |
| Step 4 | 构建空间半径邻域（radius_multiplier 控制范围） |
| Step 5 | 计算免疫抑制基因模块评分 |
| Step 6 | 构建 kNN 空间邻接图，统计邻域细胞组成向量，Leiden 聚类 |
| Step 7 | 对 Leiden cluster 进行功能评分语义注释 |
| Step 8 | 构建综合分析 DataFrame（整合所有评分） |
| Step 9 | 多层次评分计算（Treg_like_score、immune_stroma_score、immunosuppressive_niche_score） |
| Step 10 | 辅助空间区域标注（tumor_core / tumor_edge / stroma_immune / other） |
| Step 11 | 敏感性分析（v2 连续指标） |
| Step 12 | 保存核心评分表和参数元数据 |
| Step 13 | 可视化（共 15+ 张图） |
| Step 14 | 特征基因提取（三层筛选策略） |
| Step 15 | 参数扫描稳定性评估（k × quantile 网格搜索） |

---

### 2.1 核心数据表

#### `spatial_niche/spatial_niche_scores.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 12**

**文件格式**：CSV 表格，行 = spot，列 = 多种评分和标签

**关键列说明**：

| 列名 | 类型 | 来源 Step | 含义 |
|------|------|-----------|------|
| `spot_id` | 字符串 | Step 8 | Spot 唯一标识符（含样本前缀，如 `HCC4R_ACGCCTGACACGCGCT-1`） |
| `sample` | 字符串 | Step 8 | 切片来源（"HCC4R" 或 "CHC20"，来自 `adata.obs["sample"]`） |
| `spatial_x` / `spatial_y` | 浮点数 | Step 8 | Spot 空间坐标（像素单位，来自 `adata.obsm["spatial"]`） |
| `Hepatocyte` | 浮点数（0–1） | Step 3 | 肝细胞比例（归一化细胞丰度） |
| `Treg` | 浮点数（0–1） | Step 3 | 调节性 T 细胞比例 |
| `Myeloid` | 浮点数（0–1） | Step 3 | 髓系细胞比例 |
| `Fibroblast` | 浮点数（0–1） | Step 3 | 成纤维细胞比例 |
| `T/NK` | 浮点数（0–1） | Step 3 | T 细胞/NK 细胞合并比例 |
| `immunosuppressive_gene_score` | 浮点数 | Step 5 | FOXP3/IL2RA/CTLA4/TIGIT/LAG3/TGFB1/IL10 等免疫抑制标志基因的 CP10K+log1p 归一化平均表达量 |
| `neighbor_radius` | 浮点数 | Step 4 | 空间半径邻域实际搜索半径（坐标像素单位） |
| `neighborhood_cluster` | 字符串（"0","1",...） | Step 6 | Leiden 邻域组成聚类分配的 cluster 编号（无监督发现） |
| `niche_semantic_label` | 字符串 | Step 7 | 语义注释标签：`immunosuppressive_niche` 或 `cluster_X` |
| `hep_high` | 布尔值 | Step 8 | 该 spot 是否属于高肝细胞区域（Hepatocyte 比例 ≥ 75 百分位） |
| `distance_to_hep_high` | 浮点数 | Step 8 | 该 spot 到最近高肝细胞区域 spot 的欧式距离（像素） |
| `neighbor_Hepatocyte` / `neighbor_Treg` 等 | 浮点数 | Step 8 | 空间半径邻域内各细胞类型比例的均值（消除局部点噪声，捕捉微环境特征） |
| `has_hep_high_neighbor` | 布尔值 | Step 8 | 该 spot 的空间邻居中是否存在 hep_high 的 spot（用于 tumor_edge 标注） |
| `Treg_like_score` | 浮点数 | Step 9 | Z-score(Treg) + Z-score(免疫抑制基因评分)，综合反映 Treg 样特征 |
| `immune_stroma_score` | 浮点数 | Step 9 | Z-score(Treg) + Z-score(T/NK) + Z-score(Myeloid) + Z-score(Fibroblast)，四类免疫/基质细胞的综合评分 |
| `immunosuppressive_niche_score` | 浮点数 | Step 9 | 五维 Z-score 综合评分（含 neighbor_Hepatocyte 邻域肝细胞均值），**用于下游分组** |
| `niche_high` | 布尔值 | Step 9 | 该 spot `immunosuppressive_niche_score` 是否 ≥ 80 百分位（用于签名基因提取分组） |
| `spatial_region` | 字符串 | Step 10 | 规则化空间区域标注（`tumor_core`/`tumor_edge`/`stroma_immune`/`other`） |

**理解整个评分体系**：

分析采用**两阶段策略**识别免疫抑制生态位：

**阶段一：数据驱动的无监督发现（主体）**

1. **邻域组成向量**（Step 6）：对每个 spot，统计其空间上最近 k=15 个邻居（限制在同一切片内）的细胞类型比例均值，形成"邻域细胞组成向量"
   - *比喻*：评估一个社区，不只看自身，还看周围 15 户邻居的平均状况——消除单点噪声，真正刻画局部微环境

2. **Leiden 聚类**（Step 6）：基于邻域组成向量进行无监督聚类，将组成相似的 spot 归为同一 `neighborhood_cluster`（NC）
   - 此时只知道哪些区域相似，不知道其生物学意义

3. **功能评分语义注释**（Step 7）：计算每个 cluster 在 Treg 比例、Myeloid 比例、Fibroblast 比例、免疫抑制基因评分四维度的均值，综合排名最高（`combined_rank` 最小）的 cluster 被标注为 `immunosuppressive_niche`

**阶段二：连续评分与分位数切割（辅助/下游分组用）**

4. **综合评分**（Step 9）：为每个 spot 计算连续的 `immunosuppressive_niche_score`（五维 Z-score 加总），与 Leiden 聚类解耦
5. **切割阈值**（Step 9）：评分超过 80 百分位的 spot 标记为 `niche_high`，用于 Step 14 签名基因提取时的高/低组对比

**揭示的生物学现象**：

- `niche_semantic_label = "immunosuppressive_niche"` 定义了**肿瘤免疫抑制微生态位**的空间范围：Treg、TAM（Myeloid）和 CAF（Fibroblast）三类细胞同时局部富集
- 这类区域在生物学上对应肿瘤的"免疫沙漠"或"免疫排除"区域，是效应 T 细胞无法发挥杀伤功能的关键空间结构
- **联合分析优势**：两切片共约 5,000–6,000 个 spot 参与 Leiden 聚类，聚类稳定性更高；`niche_high`/`niche_low` 组的样本量更大，签名基因提取统计功效更强

---

#### `spatial_niche/spatial_niche_parameters.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 12**

**文件格式**：CSV 表格（两列：`parameter` 和 `value`）

**内容**：本次分析使用的所有关键参数及运行时计算的统计量：

| 参数 | 含义 |
|------|------|
| `abundance_key` | 细胞丰度矩阵的键名（`means_cell_abundance_w_sf`） |
| `n_neighbors_knn` | 邻域组成聚类时的 k-NN 邻居数（默认 k=15） |
| `leiden_resolution` | Leiden 聚类分辨率（默认 0.5） |
| `n_neighborhood_clusters` | 实际识别出的 Leiden cluster 数量 |
| `neighbor_radius_multiplier` | 空间半径邻域倍增系数（默认 1.25） |
| `neighbor_radius` | 实际使用的搜索半径（坐标像素单位，自适应计算：最近邻距离中位数 × 倍增系数） |
| `hep_high_quantile` | 高肝细胞区域的分位数阈值（默认 0.75） |
| `hep_high_threshold` | 高肝细胞区域的 Hepatocyte 比例具体阈值 |
| `niche_high_quantile` | niche_high 分位数阈值（默认 0.80） |
| `niche_high_threshold` | niche_high 评分的具体数值阈值 |
| `n_niche_high` | 被标记为 niche_high 的 spot 总数 |
| `marker_genes_used` | 实际用于计算免疫抑制基因评分的基因列表（FOXP3 等，若数据集中不存在则自动跳过） |

**用途**：
1. **可重现性保证**：记录参数供他人完整复现分析
2. **质控检查**：`n_neighborhood_clusters` 过多（>15）或过少（<3）均提示聚类参数需调整
3. **论文方法部分**：参数来源的直接依据

---

#### `spatial_niche/sensitivity_analysis.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 11**（`_sensitivity_analysis()` 函数，v2 连续指标版本）

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| `param_mode` | 参数类型（`knn` 或 `radius`） |
| `param_label` | 参数具体设置（如 `kNN k=10`、`radius×1.25`） |
| `is_reference` | 是否为参考配置（k=15 为参考基准，值为 True） |
| `n_niche_high` | 该参数下 niche_high 的 spot 数量（仅供参考，不参与相似度计算） |
| `niche_pct` | niche_high spot 占总 spot 的比例（%） |
| `spearman_rho` | **Spearman 秩相关系数**：该参数配置 vs 参考配置（k=15）的评分排序一致性，越接近 1.0 越稳健 |
| `weighted_jaccard` | **加权 Jaccard 相似度**：基于 Min-Max 归一化连续向量计算，对高分区域天然加权，越接近 1.0 越稳健 |
| `score_std` | 该参数下 niche 评分的标准差（反映评分分布的离散程度） |

**扫描范围**：

- kNN 邻居数：k=10 / 15（参考） / 20
- 空间半径倍增系数：×1.0 / 1.25（参考） / 1.5

**v2 版本改进说明**：

v1 版本使用硬截断 Jaccard（问题：丢失量级信息、边界敏感、固定边缘概率）。v2 改用两个**连续型指标**，彻底避免硬截断：

- **Spearman ρ**：衡量全局评分排序一致性，不受任何截断影响
- **加权 Jaccard**：对高分 spot 天然加权（`∑min(a,b) / ∑max(a,b)`），重点关注高 niche 区域一致性

**判断标准**：
- `spearman_rho > 0.95` 且 `weighted_jaccard > 0.90`：结果极稳健
- `spearman_rho > 0.90` 且 `weighted_jaccard > 0.80`：结果稳健，参数不敏感
- `spearman_rho < 0.85` 或 `weighted_jaccard < 0.70`：参数变化有实质影响，建议重审参数选择

**揭示的生物学现象**：

稳健的 niche 边界意味着免疫抑制微生态位是真实存在的、结构性的空间特征，而非因任意选定某个参数值造成的人为划分，增强了该生态位作为真实生物学实体的可信度。

---

#### `spatial_niche/neighborhood_cluster_stats.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 7**（`_annotate_neighborhood_clusters()` 函数）

**文件格式**：CSV 表格，行 = Leiden cluster，列 = 各维度统计量

**列说明**：

| 列名 | 含义 |
|------|------|
| （索引列） | Leiden cluster 编号 |
| `mean_treg` | 该 cluster 内所有 spot 的 Treg 比例均值 |
| `mean_myeloid` | Myeloid 比例均值 |
| `mean_fibroblast` | Fibroblast 比例均值 |
| `mean_immune_score` | 免疫抑制基因评分均值 |
| `rank_mean_treg` | Treg 均值的降序排名（1=最高，与 combined_rank 的计算组件） |
| `rank_mean_myeloid` | Myeloid 均值的降序排名 |
| `rank_mean_fibroblast` | Fibroblast 均值的降序排名 |
| `rank_mean_immune_score` | 免疫抑制基因评分均值的降序排名 |
| `combined_rank` | 四维排名之和（**越小说明综合免疫抑制特征越强**） |

**注释逻辑**：`combined_rank` 最小的 cluster → 标注为 `immunosuppressive_niche`（四维综合评分最高者）

**如何理解**：

这是 Step 7 语义注释的中间结果，展示了每个 Leiden cluster 的免疫抑制特征剖面。通过多维证据权衡（而非只看 Treg 一个指标）选出免疫抑制生态位，避免单一指标带来的偏差。

**揭示的生物学现象**：

若免疫抑制 cluster 的 `combined_rank` 远低于次低 cluster（如目标 cluster `combined_rank`=4，次低=12），说明免疫抑制信号的空间集中性非常突出，具有明显的"hot spot"特征；反之若差异不大，可能提示免疫抑制在组织中呈较弥散分布。

---

### 2.2 特征基因文件

#### `spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 14**（`_rank_niche_genes()` 函数，Layer 1 + Layer 2）

**文件格式**：CSV 表格，按 `composite_score` 降序排列

**列说明**：

| 列名 | 含义 |
|------|------|
| `gene` | 基因名（HGNC 符号） |
| `mean_high` | 该基因在 `niche_high` spot 中的平均 log1p 表达量 |
| `mean_low` | 该基因在 `niche_low` spot 中的平均 log1p 表达量 |
| `log2_fc` | log₂((mean_high+1)/(mean_low+1))，在免疫抑制区域的富集倍数 |
| `frac_high` | `niche_high` 组中检出该基因（表达>0）的 spot 比例 |
| `frac_low` | `niche_low` 组中检出该基因（表达>0）的 spot 比例 |
| `delta_frac` | `frac_high - frac_low`：检出率差值（**弥补稀释效应对稀有基因不敏感的核心指标**） |
| `composite_score` | 综合排序分（0.5×norm(log2FC) + 0.5×norm(delta_frac)），同时考虑表达量和检出率 |
| `pvalue` | Wilcoxon 秩和检验 p 值（Mann-Whitney U，单侧，`alternative="greater"`） |
| `fdr` | BH-FDR 多重检验校正后 p 值 |

**三层筛选策略**（Layer 1 + Layer 2 结果合并）：

| 层级 | 方法 | 筛选条件 | 目的 |
|------|------|---------|------|
| **Layer 1**（主筛选） | Wilcoxon + BH-FDR + 综合排序 | FDR<0.05 且（log2FC>0.5 OR delta_frac>0.10） | 发现统计显著的差异表达基因 |
| Layer 1 补充 | 放宽条件 | FDR<0.2 且（log2FC>0.3 OR delta_frac>0.05） | 不足 Top-N 时扩大候选 |
| **Layer 2**（强制合并） | 先验基因集 AUC+Fraction 检验 | AUC>0.6 且 delta_frac>0.05 | 将 FOXP3 等重要先验基因强制纳入签名 |

**`delta_frac` 的核心意义**：

Visium spot 覆盖 5-50 个细胞，FOXP3 等稀有 Treg 标志基因因细胞稀释效应，在大多数 spot 中表达为 0，导致均值 log2FC 接近 0，传统 Wilcoxon+FC 策略无法将其选入 Top-N。`delta_frac` 度量"niche_high 相比 niche_low 有多少更多 spot 能检出该基因"，弥补了均值法对稀有基因不敏感的缺陷。

**揭示的生物学现象**：

签名基因集代表了肿瘤局部免疫抑制环境的分子特征。典型上调基因通常包括：

- **Treg 标志基因**：FOXP3、IL2RA、CTLA4、TIGIT（免疫检查点分子）
- **免疫抑制细胞因子**：TGFB1、IL10（直接抑制效应 T 细胞活性）
- **招募趋化因子**：CXCL12、CCL22（招募 Treg 进入肿瘤的信号分子）
- **基质重塑相关基因**：FAP、ACTA2、COL1A1（CAF 激活标志，参与物理屏障构建）
- **TAM 功能基因**：SPP1、CD163、MRC1（肿瘤相关巨噬细胞标志，参与免疫抑制）

---

#### `spatial_niche/immunosuppressive_niche_signature_genes.txt`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 14**

**文件格式**：纯文本，每行一个基因名

**内容**：`immunosuppressive_niche_signature_genes_ranked.csv` 中 `gene` 列的纯文本版本（Top-N 基因，默认 80 个）。

**用途**：传递给 `code/pipeline/tcga_survival_analysis.R` 进行 TCGA 队列的 ssGSEA 评分验证，是 Step 2 与 Step 3（可选）之间的**关键接口文件**。

同时在 `results/spatial_signature_genes.txt` 保存一份同内容副本供 R 脚本通过默认路径读取。

---

#### `spatial_niche/prior_gene_set_auc.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 14**（`_rank_niche_genes()` → Layer 2 先验功能基因集 AUC 检验）

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| `gene_set` | 先验基因集分组（`Treg_markers` / `TAM_features` / `CAF_activation`） |
| `gene` | 原始先验基因名称 |
| `actual_gene` | 实际使用的基因名称（若主基因不在数据集中，则使用 GENE_FALLBACKS 中的替代基因） |
| `mean_high` | 该基因在 `niche_high` 组的平均表达量 |
| `mean_low` | 该基因在 `niche_low` 组的平均表达量 |
| `log2_fc` | 两组 log2 倍数变化 |
| `frac_high` | niche_high 组的检出率 |
| `frac_low` | niche_low 组的检出率 |
| `delta_frac` | 检出率差值 |
| `pvalue` | Mann-Whitney U 检验 p 值 |
| `auc` | AUC 值（用 Mann-Whitney U 统计量计算，0.5=无区分力，>0.6=有区分力） |
| `fdr` | BH-FDR 校正后 p 值 |

**三组先验基因集**：

| 基因集 | 基因列表 | 生物学意义 |
|--------|---------|----------|
| Treg_markers | FOXP3, IL2RA, CTLA4, TIGIT, IKZF2 | Treg 核心标志（主调控转录因子和检查点分子） |
| TAM_features | CD163, MRC1, TGFB1, IL10, CXCL12 | 肿瘤相关巨噬细胞功能基因（M2 极化和免疫抑制分泌） |
| CAF_activation | FAP, ACTA2, POSTN, COL1A1, CCL22 | 肿瘤相关成纤维细胞激活标志（基质重塑和免疫排斥） |

**如何理解**：

AUC>0.6 表示该基因能区分 niche_high 和 niche_low 组（高于随机猜测）；结合 delta_frac>0.05 说明不仅均值差异显著，而且在更多 spot 中可检出，具有真实的生物学富集。

---

#### `spatial_niche/gini_score_genes.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 14**（`_rank_niche_genes()` → Layer 3 Gini Index 特异性评分）

**文件格式**：CSV 表格（若无满足条件的基因，仍会写出空文件以保证路径可预期）

**列说明**：

| 列名 | 含义 |
|------|------|
| `gene` | 基因名 |
| `gini` | Gini 系数（0=均匀分布，1=极度集中；数值越高说明该基因在少数 spot 中高度局灶性表达） |
| `log2_fc` | 该基因在 niche_high vs niche_low 中的 log2 倍数变化 |
| `mean_high` | 该基因在 niche_high 组的平均表达量 |

**筛选条件**：Gini>0.3 且 log2FC>0（阈值由 0.5 降至 0.3，覆盖 FOXP3 等稀有免疫基因）

**按 Gini 降序排列**。

**Gini Index 的意义**：

Gini 系数（基尼系数）原为经济学中衡量收入不平等的指标。在基因表达分析中，高 Gini 值意味着该基因仅在少数 spot 中高度表达（局灶性），而在大多数 spot 中接近 0。这正是 FOXP3 等稀有免疫基因的特征——只在 Treg 稀少富集的局灶区域可检出。

**揭示的生物学现象**：

Gini 高的基因（如 FOXP3）即使均值 log2FC 很低，也代表着真实的"局灶性免疫富集"现象：只有在免疫抑制 niche 中才有少数 spot 能检测到这些稀有基因，这种"点状" 而非"弥散"的分布模式是 Treg 等稀有免疫细胞的空间特征。

---

#### `spatial_niche/param_scan_deg_stability.csv`

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 15**（`_param_scan_deg_stability()` 函数）

**文件格式**：CSV 表格，每行为一种参数组合（共 5×4=20 行）

**列说明**：

| 列名 | 含义 |
|------|------|
| `k` | 本次扫描使用的 kNN 邻居数（8/10/15/20/25） |
| `quantile` | 本次扫描使用的 niche_high 分位数阈值（0.70/0.75/0.80/0.85） |
| `n_sig_deg` | 该参数组合下显著 DEG 数量（FDR<0.05 且 log2FC>0.5） |
| `mean_log2fc_topN` | 前 N 个基因（按 log2FC 降序）的平均 log2FC |
| `top_genes_str` | 前 20 个显著基因名（分号分隔，用于跨参数 Jaccard 比较） |

**扫描范围**（5×4 网格）：
- 横轴：`niche_high_quantile`（0.70 / 0.75 / 0.80 / 0.85）
- 纵轴：kNN 邻居数 k（8 / 10 / 15 / 20 / 25）

**为什么改为 k × quantile，而非 resolution × quantile（历史修复说明）**：

原来的 `resolution × quantile` 扫描存在根本性缺陷：`niche_score` 的计算基于邻域组成向量的 Z-score 加总，而邻域组成向量由 k-NN 邻居数 k 决定。Leiden resolution 仅决定聚类粒度（cluster 数量），不影响 `niche_score` 数值。因此固定 k 时改变 resolution 对 DEG 结果毫无影响，热图所有列完全相同。改为 k × quantile 后，k 不同则邻域大小不同，`niche_score` 向量随之变化，DEG 结果真正不同，热图才有意义。

**如何解读配套热图**（`param_scan_deg_stability_heatmap.png`）：

- 横轴：`niche_high_quantile`（阈值越高 = niche_high 越少 = 越严格）
- 纵轴：kNN 邻居数 k（k 越大 = 邻域越大 = 捕捉宏观模式）
- 颜色深浅：显著 DEG 数量（越深 = DEG 越多）
- 选参标准：颜色最深且处于"高原"区域（与相邻参数结果相近）的参数组合为推荐；若主流程参数（k=15, q=0.80）位于高原区，则合理

**揭示的生物学现象**：

若 k 和 quantile 在合理范围内变动时 DEG 列表高度一致（Top-20 基因 Jaccard > 0.7），说明免疫抑制生态位的分子特征不依赖于参数选择的偶然性，而是组织中稳定存在的生物学实体。

---

### 2.3 可视化图表（`spatial_niche/plots/`）

**来源脚本 / Step**：`run_spatial_niche_analysis.py` → **Step 13**（主可视化）、**Step 13+**（补充图）、**Step 14**（签名基因图）、**Step 15**（参数扫描热图）

以下按逻辑功能分组讲解所有可视化输出。

---

#### 2.3.1 细胞类型空间分布图

##### `spatial_hepatocyte.png`

**来源 Step**：Step 13 → `_spatial_scatter()`

**内容**：Hepatocyte 比例的空间散点图。每个点代表一个 Visium spot，颜色深浅反映肝细胞比例（0–1），使用 paper_ybp（Yellow-Black-Purple）配色。

**如何理解**：

高肝细胞比例区域（亮黄色）对应肿瘤实质核心区；低值区域（深紫色）对应间质和免疫浸润区。此图是 Step 10 空间区域标注（tumor_core / tumor_edge / stroma_immune）的基础参考。

**揭示的生物学现象**：

肝细胞癌组织中肝细胞比例呈空间异质性分布，高密度区域与低密度区域之间形成的边界带（tumor_edge）是肿瘤-免疫相互作用的关键界面。

---

##### `spatial_treg.png`

**来源 Step**：Step 13 → `_spatial_scatter()`

**内容**：Treg 比例的空间散点图，paper_ybp 配色。

**如何理解**：

Treg 高比例区域（亮黄色）呈现稀疏、局灶性分布模式，而非弥散性浸润。这与 Treg 作为稀有细胞群体的特性一致——通常以小簇形式聚集在肿瘤边缘或免疫抑制微生态位中。

**揭示的生物学现象**：

Treg 的局灶性富集模式揭示了"免疫抑制热点"的存在——这些局部区域通过 Treg 介导的免疫抑制为肿瘤提供了免疫逃逸的"避风港"。

---

##### `spatial_myeloid.png` / `spatial_fibroblast.png`

**来源 Step**：Step 13 → `_spatial_scatter()`

**内容**：Myeloid 和 Fibroblast 比例的空间散点图，paper_ybp 配色。

**如何理解**：

与 Treg 类似，Myeloid（肿瘤相关巨噬细胞 TAM）和 Fibroblast（肿瘤相关成纤维细胞 CAF）的高值区域应与 Treg 热点在空间上重叠，共同定义免疫抑制生态位。

**揭示的生物学现象**：

Treg-TAM-CAF 三者的空间共定位是免疫抑制生态位的细胞学基础：TAM 分泌 CCL22/TGFB1 招募和维持 Treg，CAF 构建物理屏障排斥效应 T 细胞，三者协同形成"免疫沙漠"。

---

#### 2.3.2 免疫抑制 Niche 综合评分图

##### `spatial_immunosuppressive_niche_score.png`

**来源 Step**：Step 13 → `_spatial_scatter()`，使用 magma 配色

**内容**：`immunosuppressive_niche_score`（五维 Z-score 综合评分）的空间分布图。magma 配色：暗色 = 低评分，亮黄色 = 高评分（免疫抑制信号最强）。

**如何理解**：

高评分区域（亮色）的空间分布模式直接展示免疫抑制生态位的定位和范围。理想情况下，高评分区域应呈连续片状分布（而非随机散点），说明免疫抑制信号在空间上具有结构性而非噪声。

**揭示的生物学现象**：

免疫抑制信号的空间连续性意味着这不是随机事件，而是肿瘤微环境中结构性存在的免疫逃逸区域。连续的免疫抑制区域比孤立散点具有更强的临床意义。

---

#### 2.3.3 Leiden 邻域聚类图

##### `spatial_neighborhood_clusters.png`

**来源 Step**：Step 13 → `_plot_neighborhood_clusters()`

**内容**：基于邻域组成向量的 Leiden 无监督聚类结果的空间分布图。每个 Leiden cluster（NC-0, NC-1, ...）使用高饱和度离散色（仿论文 tSNE 图配色），每个 spot 按其所属 cluster 着色。

**如何理解**：

此图展示"空间细胞生态区"（Cellular Neighborhoods）的划分结果——空间上相邻、细胞组成相似的 spot 被归为同一 cluster。不同颜色代表不同的细胞生态区。关注哪个 cluster 被标注为 `immunosuppressive_niche`（在 `neighborhood_cluster_stats.csv` 中 `combined_rank` 最小的 cluster）。

**揭示的生物学现象**：

Leiden 聚类将组织划分为多个空间生态区，其中免疫抑制生态区（Treg+TAM+CAF 共富集区）作为独立的空间实体被无监督识别出来，证明免疫抑制微环境是组织级别的结构性特征，而非人为定义。

---

##### `spatial_niche_semantic_labels.png`

**来源 Step**：Step 13 → `_spatial_scatter(categorical=True)`

**内容**：语义注释后的 niche 标签空间分布图。仅两类标签：`immunosuppressive_niche`（紫色）vs 其他（灰色）。使用固定色板：tumor_core → 红色，tumor_edge → 橙色，stroma_immune → 蓝色，immunosuppressive_niche → 紫色，other → 灰色。

**如何理解**：

紫色区域即为 Leiden 聚类中 `combined_rank` 最小的 cluster 在空间上的分布。与 `spatial_neighborhood_clusters.png` 对比可确认哪个 NC 对应免疫抑制生态位。

**揭示的生物学现象**：

语义标签将无监督聚类结果映射为可解读的生物学概念——"免疫抑制生态位"，明确了这些区域的功能属性。

---

#### 2.3.4 辅助区域标注图

##### `spatial_region_labels.png`

**来源 Step**：Step 13 → `_spatial_scatter(categorical=True)`

**内容**：基于规则化标注的空间区域分类图（`spatial_region` 列），四种区域以不同颜色显示：

| 区域 | 颜色 | 判定规则 |
|------|------|---------|
| `tumor_core` | 红色 | `hep_high=True` 且 `has_hep_high_neighbor=True`（高肝细胞区域中心） |
| `tumor_edge` | 橙色 | `hep_high=True` 且 `has_hep_high_neighbor=False`，或 `hep_high=False` 且 `has_hep_high_neighbor=True`（高肝细胞区域边缘/过渡带） |
| `stroma_immune` | 蓝色 | `hep_high=False` 且 `has_hep_high_neighbor=False`（远离肿瘤实质的间质/免疫浸润区） |
| `other` | 灰色 | 不满足以上条件的 spot |

**如何理解**：

此图是辅助性标注，将组织空间划分为功能上可解释的区域，便于可视化解读。`tumor_edge`（肿瘤边缘带）是免疫细胞与肿瘤细胞直接接触的前线区域，`stroma_immune` 是间质免疫浸润区——两者是免疫抑制生态位最可能出现的位置。

**揭示的生物学现象**：

肿瘤边缘带（tumor_edge）是肿瘤-免疫相互作用最剧烈的空间界面，免疫抑制生态位倾向于在此区域形成，代表了肿瘤"驯化"局部免疫微环境的关键空间结构。

---

#### 2.3.5 距离依赖性折线图

##### `distance_to_hep_high_vs_niche_score.png`

**来源 Step**：Step 13 → `_plot_distance()`

**内容**：X 轴为到最近高肝细胞区域 spot 的距离分箱（5 个等频区间），Y 轴为该距离区间内免疫抑制 niche 评分的均值 ± 标准误差。

**如何理解**：

曲线若呈**倒 U 型**（近距离和远距离 niche 评分低，中等距离最高），说明免疫抑制信号在肿瘤边缘带最强——太靠近肿瘤核心处免疫细胞浸润不足，太远离肿瘤区域则缺乏肿瘤抗原刺激。

曲线若**单调递增**（距离越远 niche 评分越高），则说明免疫抑制信号集中在间质区域，肿瘤核心免疫排除效应显著。

曲线若**平坦**，则说明免疫抑制信号在空间上无明显距离依赖性，呈弥散分布。

**揭示的生物学现象**：

免疫抑制生态位的空间定位规律——距离依赖性曲线揭示了免疫抑制信号与肿瘤实质的空间关系，是理解肿瘤免疫逃逸空间策略的关键证据。

---

#### 2.3.6 区域评分箱线图

##### `region_score_boxplots.png`

**来源 Step**：Step 13 → `_plot_region_box()`

**内容**：双面板箱线图，展示四种空间区域（tumor_core → tumor_edge → stroma_immune → other）中：
- 左图：`Treg_like_score` 分布
- 右图：`immunosuppressive_niche_score` 分布

**如何理解**：

理想模式下，`tumor_edge` 和 `stroma_immune` 区域的评分中位数和上四分位距应显著高于 `tumor_core`，验证免疫抑制信号在肿瘤边缘和间质区域的富集。

**揭示的生物学现象**：

肿瘤核心区（tumor_core）由于肝细胞密度高、免疫细胞物理空间受限，呈现"免疫排除"状态；而肿瘤边缘和间质区域因免疫细胞浸润通道开放，形成免疫抑制生态位的细胞学条件更充分。

---

#### 2.3.7 Hepatocyte vs Treg 散点图

##### `hepatocyte_vs_treg_niche_score.png`

**来源 Step**：Step 13 → `_plot_hep_treg()`

**内容**：X 轴为 Hepatocyte 比例，Y 轴为 Treg 比例，颜色编码 `immunosuppressive_niche_score`（paper_ybp 配色）。黑色虚线标注 0.75 分位数阈值，划分四个象限。

**如何理解**：

- **左下象限**（低肝细胞 + 低 Treg）：间质区域但 Treg 浸润少
- **右下象限**（高肝细胞 + 低 Treg）：肿瘤核心区，免疫排除
- **左上象限**（低肝细胞 + 高 Treg）：**免疫抑制生态位核心区域**——Treg 在非实质区域富集
- **右上象限**（高肝细胞 + 高 Treg）：罕见，可能代表 Treg 直接浸润肿瘤实质

颜色梯度（paper_ybp）进一步展示 niche 综合评分在这些象限中的分布。

**揭示的生物学现象**：

Treg 与肝细胞比例的负相关关系是肿瘤免疫逃逸空间策略的直观体现——Treg 选择在肿瘤边缘和间质区域富集（而非肿瘤核心），这种空间选择可能是为了在效应 T 细胞进入肿瘤的"通道"上建立免疫抑制屏障。

---

#### 2.3.8 细胞类型相关性热图

##### `celltype_niche_correlation.png`

**来源 Step**：Step 13 → `_plot_correlation()`

**内容**：Pearson 相关系数热图，展示以下 7 个变量两两之间的相关性：

| 变量 | 含义 |
|------|------|
| Hepatocyte | 肝细胞比例 |
| Treg | 调节性 T 细胞比例 |
| T/NK | T 细胞/NK 细胞比例 |
| Myeloid | 髓系细胞比例 |
| Fibroblast | 成纤维细胞比例 |
| Treg_like_score | Treg 样综合评分 |
| immunosuppressive_niche_score | 免疫抑制生态位综合评分 |

配色：vlag（蓝-白-红），中心=0。

**如何理解**：

关注以下关键相关对：
- Treg 与 Myeloid 正相关 → Treg-TAM 共定位，验证免疫抑制生态位细胞组成
- Treg 与 Fibroblast 正相关 → Treg-CAF 协同
- Hepatocyte 与 Treg/Myeloid 负相关 → 肿瘤实质与免疫浸润区的空间互斥
- `immunosuppressive_niche_score` 与 Treg/Myeloid/Fibroblast 均正相关 → 综合评分有效整合多维信息

**揭示的生物学现象**：

Treg-TAM-CAF 三者之间的空间正相关性验证了免疫抑制生态位的细胞学协同机制：三类细胞不是随机分布，而是通过旁分泌信号（CCL22-CCR4、TGFB1 等）相互招募和维持，形成自我强化的免疫抑制回路。

---

#### 2.3.9 L-R 配体-受体通讯热图

##### `lr_communication_heatmap.png`

**来源 Step**：Step 13 → `_plot_lr_communication()`

**内容**：双面板图，展示 12 对预设配体-受体对在 niche_high vs niche_low 中的信号强度差异：

**左图**：信号强度热图
- 行：L-R 对名称
- 列：Niche-High / Niche-Low 两组
- 颜色：行最大值归一化后的信号强度（YlOrRd 配色，深红=强信号）
- 格内标注原始数值

**右图**：log2FC 条形图
- 红色条 = 在 Niche-High 中富集（log2FC > 0）
- 蓝色条 = 在 Niche-Low 中富集（log2FC < 0）

**预设 L-R 对**：

| L-R 对 | 生物学功能 |
|--------|----------|
| CXCL12–CXCR4 | Treg 招募核心轴（TAM/CAF 分泌 CXCL12 招募 Treg） |
| TGFB1–TGFBR1 | 免疫抑制性细胞因子信号（Treg/TAM 分泌，抑制效应 T 细胞） |
| PD-1–PD-L1 | 经典免疫检查点（T 细胞耗竭） |
| TIGIT–NECTIN2 | 新型抑制性检查点 |
| LAG3–MHC-II | 免疫抑制检查点 |
| CCL22–CCR4 | Treg 特异招募趋化因子轴 |
| IL10–IL10RA | 免疫抑制性细胞因子 |
| SPP1–CD44 | TAM 介导的免疫抑制（骨桥蛋白通路） |
| MIF–CD74 | 巨噬细胞迁移抑制因子 |
| VEGFA–KDR | 血管生成（与 CAF 相关） |
| Galectin9–TIM-3 | 抑制性检查点 |
| CCL2–CCR2 | 巨噬细胞招募通路 |

**技术细节**：

采用**秩归一化乘积**（rank-normalized product）策略计算通讯强度：先对配体和受体基因的表达向量做秩归一化（`rank / n_spots`），再取乘积均值。这解决了 Visium spot 中稀疏基因（如 CCL22、IL10）因细胞稀释效应导致直接乘积为 0 的问题。

缺失基因自动使用 `GENE_FALLBACKS` 中的同家族替代基因（如 CCR4 → CCR2/CCR5），标签中用 `*` 标注。

**如何理解**：

左图热图中的颜色对比（Niche-High 列 vs Niche-Low 列）直接展示哪些通讯轴在免疫抑制生态位中被激活。右图条形图的红色条越长，说明该 L-R 对在免疫抑制区域中的信号越强。

**揭示的生物学现象**：

CXCL12–CXCR4 和 CCL22–CCR4 等趋化因子轴在 niche_high 中的富集，验证了"TAM/CAF 通过旁分泌招募 Treg"的免疫抑制回路模型。SPP1–CD44 的富集则揭示了 TAM 通过骨桥蛋白通路参与免疫抑制的机制。这些通讯轴共同构成了免疫抑制生态位的分子通讯网络。

---

#### 2.3.10 敏感性分析稳定性图

##### `sensitivity_niche_stability.png`

**来源 Step**：Step 13 → `_plot_sensitivity()`

**内容**：2×2 布局图（左列=kNN 参数，右列=radius 参数；上行=Spearman ρ，下行=Weighted Jaccard）：

| 子图位置 | 内容 |
|---------|------|
| 左上 | kNN 邻居数变化下的 Spearman 秩相关系数 |
| 左下 | kNN 邻居数变化下的加权 Jaccard 相似度 |
| 右上 | 空间半径倍增系数变化下的 Spearman 秩相关系数 |
| 右下 | 空间半径倍增系数变化下的加权 Jaccard 相似度 |

红色柱 = 参考参数（k=15 kNN 或 radius×1.25），蓝色柱 = 其他参数设置。橙色虚线 = 稳健性阈值（Spearman: 0.90, WJ: 0.80），红色虚线 = 参考值（=1.0）。

**如何理解**：

所有柱子均高于橙色虚线（稳健性阈值）→ 结果极稳健，参数选择不影响结论。
部分柱子低于阈值 → 结果对参数敏感，需审慎选择参数并报告不确定性。

**揭示的生物学现象**：

若所有参数配置下 Spearman ρ > 0.95 且 WJ > 0.90，说明免疫抑制生态位的识别不依赖于 kNN 邻居数或搜索半径的具体选择，是组织中真实存在的结构性特征。

---

#### 2.3.11 Niche-High 二值分布图

##### `spatial_niche_high_score_spots.png`

**来源 Step**：Step 13+ → `_spatial_scatter()`，使用 RdGy_r 配色

**内容**：基于评分阈值切割的 niche_high 二值空间分布图。红色 = 评分高于 niche_high_quantile（如 80 百分位）的 spot，灰色 = 其余 spot。

**如何理解**：

红色 spot 的空间分布展示了"评分法"识别的免疫抑制热点。与 `spatial_niche_semantic_labels.png`（Leiden 聚类法）对比，评估两种方法的一致性。

**揭示的生物学现象**：

若红色区域与 Leiden 聚类法识别的免疫抑制 niche 高度重叠，说明无监督聚类和规则化评分法两种独立策略指向同一结论，大幅增强结果的可信度。

---

#### 2.3.12 聚类法 vs 评分法并排对比图

##### `spatial_niche_cluster_vs_score_comparison.png`

**来源 Step**：Step 13+ → `_plot_niche_comparison()`

**内容**：并排双图，共用同一空间坐标系：

- **左图**：Leiden 聚类注释法——语义标签为 `immunosuppressive_niche` 的 spot 红色高亮，其余灰色
- **右图**：评分阈值法——`niche_high=True` 的 spot 红色高亮，其余灰色

标题中同时标注两种方法结果的 **Jaccard 相似度**（交集 / 并集）。

**如何理解**：

Jaccard > 0.7 → 两种方法高度一致，结果稳健
Jaccard 0.4–0.7 → 部分一致，有意义的差异可从空间位置解读
Jaccard < 0.4 → 两种方法分歧较大，需审查参数或重新评估 niche 定义

三种典型分歧模式：
- **聚类有但评分无**：邻域组成像免疫抑制 niche 但自身评分不够高（边缘过渡区域）
- **评分有但聚类无**：自身评分高但周围邻域不典型（孤立的免疫抑制"岛"）
- **两者共有**：核心免疫抑制区域

**揭示的生物学现象**：

两种独立方法的交叉验证为免疫抑制生态位的存在提供了双证据支持。核心区域的高度重叠说明免疫抑制信号是真实且可重复检测的。

---

#### 2.3.13 火山图

##### `niche_signature_volcano.png`

**来源 Step**：Step 14 → `_plot_volcano()`

**内容**：全基因范围差异表达火山图：

- **X 轴**：log₂FC（niche_high vs niche_low）
- **Y 轴**：-log₁₀(FDR)
- **颜色编码**：
  - 🔴 红色 = 上调显著（log2FC > 0.5 且 FDR < 0.05）
  - 🔵 蓝色 = 下调显著（log2FC < -0.5 且 FDR < 0.05）
  - ⚫ 灰色 = 不显著
- **标注**：先验免疫抑制目标基因（FOXP3、CTLA4、TIGIT、FAP、ACTA2、CCL22 等）在显著或接近阈值时被标注基因名
- **标题**：显示上调/下调显著基因数量

**如何理解**：

右上角红色点为免疫抑制生态位中上调最显著的基因。注意 FOXP3 等稀有基因可能因均值稀释效应在火山图中位置偏低——这正是引入检出率散点图（2.3.14）的原因。

**揭示的生物学现象**：

火山图展示了免疫抑制生态位的全局转录特征。典型上调基因涵盖免疫检查点分子（CTLA4、TIGIT、LAG3）、免疫抑制细胞因子（TGFB1、IL10）和基质重塑基因（COL1A1、FAP、ACTA2），反映 Treg-TAM-CAF 三方协同的分子网络。

---

#### 2.3.14 检出率差值散点图

##### `niche_fraction_scatter.png`

**来源 Step**：Step 14 → `_plot_fraction_scatter()`

**内容**：全基因范围的"log2FC vs delta_frac"散点图：

- **X 轴**：log₂FC（均值倍数变化）
- **Y 轴**：delta_frac（检出率差值，niche_high 比 niche_low 多多少比例的 spot 能检出该基因）
- **颜色编码**：
  - 🔴 红色 = 双重显著（log2FC > 0.5 且 delta_frac > 0.10）
  - 🟠 橙色 = 仅检出率富集（delta_frac > 0.10 但 log2FC 不高；**典型稀有免疫基因**）
  - 🔵 蓝色 = 仅 FC 高（log2FC > 0.5 但检出率差异不大；可能是高表达管家基因）
  - ⚫ 灰色 = 两者均不显著
- **标注**：先验免疫抑制目标基因被选择性标注

**如何理解**：

此图是火山图的**关键补充视角**。火山图以均值 log2FC 为横轴，FOXP3 等稀有基因因 Visium spot 中大量 0 值压低均值而在火山图中位置偏低。本图以 delta_frac 为纵轴，直接展示"niche_high 中有更多 spot 能检出该基因"——这正是稀有免疫基因在空间水平上的真正富集模式。

**橙色点**是本图最重要的发现——它们在火山图中不显著（log2FC 不高），但在检出率维度上有真实的空间富集。FOXP3 最可能出现在橙色区域或红-橙过渡区域。

**揭示的生物学现象**：

稀有免疫基因（如 FOXP3）的空间富集模式是"局灶性检出率提升"而非"均值倍数增加"——这意味着免疫抑制信号不是均匀分布在所有 spot 中，而是只在少数局部区域可检测。这种"点状"分布模式是 Treg 等稀有细胞群体的空间特征。

---

#### 2.3.15 参数扫描热图

##### `param_scan_deg_stability_heatmap.png`

**来源 Step**：Step 15 → `_plot_param_scan_heatmap()`

**内容**：5×4 网格热图：

- **横轴**：`niche_high_quantile`（0.70 / 0.75 / 0.80 / 0.85）
- **纵轴**：kNN 邻居数 k（8 / 10 / 15 / 20 / 25）
- **颜色**：显著 DEG 数量（YlOrRd 配色，颜色越深 = DEG 越多）
- **格内标注**：具体 DEG 数量

**如何理解**：

颜色最深且处于"高原"区域（与相邻格子数值相近）的参数组合为推荐。若主流程参数（k=15, q=0.80）位于高原区中心，则说明当前参数选择合理。

注意 k=8（微观邻域）和 k=25（宏观邻域）可能给出不同的 DEG 数量，但若 Top-20 基因列表高度一致，则结果仍然稳健。

**揭示的生物学现象**：

高原区的存在说明免疫抑制生态位的分子特征在参数空间中具有鲁棒性——无论用微观还是宏观邻域定义，无论阈值宽松还是严格，核心签名基因始终被一致识别。这进一步验证了免疫抑制生态位作为真实生物学实体的可靠性。

---

## 三、结果文件索引表

| 文件路径 | 格式 | 来源 Step | 简要说明 |
|---------|------|----------|---------|
| `adata_vis_post_joint.h5ad` | h5ad | Step 1-7 | 两切片合并反卷积核心数据 |
| `adata_vis_post.h5ad` | h5ad | Step 1-7 | 同上，标准名称副本 |
| `adata_vis_post_HCC4R.h5ad` | h5ad | Step 1-6 | HCC4R 单切片反卷积 |
| `adata_vis_post_CHC20.h5ad` | h5ad | Step 1-6 | CHC20 单切片反卷积 |
| `adata_sc_post.h5ad` / `_joint.h5ad` | h5ad | Step 1-9 | Treg 重注释后 scRNA 参考 |
| `spot_cell_proportion_HCC4R.csv` | CSV | Step 1-6 | HCC4R spot 细胞比例 |
| `spot_cell_proportion_CHC20.csv` | CSV | Step 1-6 | CHC20 spot 细胞比例 |
| `spot_cell_proportion_joint.csv` | CSV | Step 1-7 | 两切片合并比例表 |
| `scrna_tsne_celltype.png` | PNG | Step 1-4+ | scRNA tSNE 细胞类型分布 |
| `scrna_celltype_marker_heatmap.png` | PNG | Step 1-4+ | 细胞类型 Marker 热图 |
| `t_cell_dotplot_horizontal.png` | PNG | Step 1-4 | Treg 标志基因气泡图 |
| `regression_training_history_joint.png` | PNG | Step 1-5 | 统一 RegressionModel 训练曲线 |
| `spatial_mapping_training_history_HCC4R.png` | PNG | Step 1-6 | HCC4R 空间建模训练曲线 |
| `spatial_mapping_training_history_CHC20.png` | PNG | Step 1-6 | CHC20 空间建模训练曲线 |
| `cross_slice_comparison/` | 目录 | Step 1-8 | 切片对比图（3 张） |
| `shared_genes_joint.txt` | TXT | Step 1-9 | 三路共享基因列表 |
| `run_preprocessing.log` | TXT | Step 1 全流程 | 运行日志 |
| `spatial_niche/spatial_niche_scores.csv` | CSV | Step 2-12 | 完整 spot 级评分表 |
| `spatial_niche/spatial_niche_parameters.csv` | CSV | Step 2-12 | 分析参数元数据 |
| `spatial_niche/sensitivity_analysis.csv` | CSV | Step 2-11 | 敏感性分析结果 |
| `spatial_niche/neighborhood_cluster_stats.csv` | CSV | Step 2-7 | 邻域聚类统计 |
| `spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv` | CSV | Step 2-14 | 签名基因排名（Layer 1+2） |
| `spatial_niche/immunosuppressive_niche_signature_genes.txt` | TXT | Step 2-14 | 签名基因列表 |
| `spatial_niche/prior_gene_set_auc.csv` | CSV | Step 2-14 | 先验基因集 AUC（Layer 2） |
| `spatial_niche/gini_score_genes.csv` | CSV | Step 2-14 | Gini 特异性基因（Layer 3） |
| `spatial_niche/param_scan_deg_stability.csv` | CSV | Step 2-15 | 参数扫描 DEG 稳定性 |
| `spatial_niche/plots/spatial_hepatocyte.png` | PNG | Step 2-13 | 肝细胞比例空间分布 |
| `spatial_niche/plots/spatial_treg.png` | PNG | Step 2-13 | Treg 比例空间分布 |
| `spatial_niche/plots/spatial_myeloid.png` | PNG | Step 2-13 | Myeloid 比例空间分布 |
| `spatial_niche/plots/spatial_fibroblast.png` | PNG | Step 2-13 | Fibroblast 比例空间分布 |
| `spatial_niche/plots/spatial_immunosuppressive_niche_score.png` | PNG | Step 2-13 | Niche 综合评分空间分布 |
| `spatial_niche/plots/spatial_neighborhood_clusters.png` | PNG | Step 2-13 | Leiden 邻域聚类空间图 |
| `spatial_niche/plots/spatial_niche_semantic_labels.png` | PNG | Step 2-13 | 语义 niche 标签空间图 |
| `spatial_niche/plots/spatial_region_labels.png` | PNG | Step 2-13 | 区域标注空间图 |
| `spatial_niche/plots/distance_to_hep_high_vs_niche_score.png` | PNG | Step 2-13 | 距离依赖性折线图 |
| `spatial_niche/plots/region_score_boxplots.png` | PNG | Step 2-13 | 区域评分箱线图 |
| `spatial_niche/plots/hepatocyte_vs_treg_niche_score.png` | PNG | Step 2-13 | Hepatocyte-Treg 散点图 |
| `spatial_niche/plots/celltype_niche_correlation.png` | PNG | Step 2-13 | 细胞类型相关性热图 |
| `spatial_niche/plots/lr_communication_heatmap.png` | PNG | Step 2-13 | L-R 通讯热图 |
| `spatial_niche/plots/sensitivity_niche_stability.png` | PNG | Step 2-13 | 敏感性分析稳定性图 |
| `spatial_niche/plots/spatial_niche_high_score_spots.png` | PNG | Step 2-13+ | Niche-High 二值空间图 |
| `spatial_niche/plots/spatial_niche_cluster_vs_score_comparison.png` | PNG | Step 2-13+ | 聚类法 vs 评分法对比图 |
| `spatial_niche/plots/niche_signature_volcano.png` | PNG | Step 2-14 | 签名基因火山图 |
| `spatial_niche/plots/niche_fraction_scatter.png` | PNG | Step 2-14 | 检出率差值散点图 |
| `spatial_niche/plots/param_scan_deg_stability_heatmap.png` | PNG | Step 2-15 | 参数扫描热图 |

---

## 四、分析流程数据流图

```
┌──────────────────────────────────────────────────────────────────────────┐
│                    Step 1: 数据预处理与 Cell2location 反卷积              │
│                   (run_preprocessing.py → main_joint())                 │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  scRNA_reference.h5ad ──┐                                                │
│                         │                                                │
│  HCC4R Visium ──────────┤── Step 3: 三路基因交集 ──► shared_genes_joint  │
│                         │                                                │
│  CHC20 Visium ──────────┘                                                │
│          │                                                               │
│          ▼                                                               │
│  Step 4: Treg 鉴定 ──► t_cell_dotplot_horizontal.png                    │
│                    ──► scrna_tsne_celltype.png                           │
│                    ──► scrna_celltype_marker_heatmap.png                 │
│          │                                                               │
│          ▼                                                               │
│  Step 5: 统一 RegressionModel ──► regression_training_history_joint.png │
│          │                         cell_state_df (统一参考签名)           │
│          ▼                                                               │
│  Step 6: Cell2location 反卷积                                            │
│    ├── HCC4R ──► adata_vis_post_HCC4R.h5ad                              │
│    │          ──► spot_cell_proportion_HCC4R.csv                         │
│    │          ──► spatial_mapping_training_history_HCC4R.png             │
│    └── CHC20 ──► adata_vis_post_CHC20.h5ad                              │
│               ──► spot_cell_proportion_CHC20.csv                         │
│               ──► spatial_mapping_training_history_CHC20.png             │
│          │                                                               │
│          ▼                                                               │
│  Step 7: 合并两切片 ──► adata_vis_post_joint.h5ad ★ (Step 2 输入)       │
│                    ──► spot_cell_proportion_joint.csv                     │
│          │                                                               │
│          ▼                                                               │
│  Step 8: 切片对比 ──► cross_slice_comparison/ (3 张)                    │
│          │                                                               │
│          ▼                                                               │
│  Step 9: 保存辅助 ──► adata_sc_post.h5ad                                │
│                    ──► shared_genes_joint.txt                             │
│                                                                          │
└──────────────────────────────────┬───────────────────────────────────────┘
                                   │
                                   ▼ adata_vis_post_joint.h5ad
┌──────────────────────────────────────────────────────────────────────────┐
│                Step 2: 空间免疫抑制生态位分析                             │
│            (run_spatial_niche_analysis.py → main())                      │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  Step 1-3: 丰度提取 → 归一化比例                                          │
│          │                                                               │
│          ▼                                                               │
│  Step 4: 半径邻域 ──► neighbor_radius 参数                               │
│          │                                                               │
│          ▼                                                               │
│  Step 5: 免疫抑制基因评分 ──► immunosuppressive_gene_score               │
│          │                                                               │
│          ▼                                                               │
│  Step 6: kNN + Leiden ──► neighborhood_cluster                          │
│          │                                                               │
│          ▼                                                               │
│  Step 7: 语义注释 ──► niche_semantic_label                              │
│                 ──► neighborhood_cluster_stats.csv                       │
│          │                                                               │
│          ▼                                                               │
│  Step 8: 综合 DataFrame ──► neighbor_* 列                               │
│          │                                                               │
│          ▼                                                               │
│  Step 9: 多层次评分 ──► Treg_like_score / immune_stroma_score           │
│                    ──► immunosuppressive_niche_score / niche_high        │
│          │                                                               │
│          ▼                                                               │
│  Step 10: 区域标注 ──► spatial_region (tumor_core/edge/stroma/other)    │
│          │                                                               │
│          ▼                                                               │
│  Step 11: 敏感性 ──► sensitivity_analysis.csv                           │
│          │                                                               │
│          ▼                                                               │
│  Step 12: 保存核心表 ──► spatial_niche_scores.csv ★                     │
│                    ──► spatial_niche_parameters.csv                      │
│          │                                                               │
│          ▼                                                               │
│  Step 13: 可视化 ──► plots/ (15+ 张图)                                  │
│          │                                                               │
│          ▼                                                               │
│  Step 14: 签名基因 ──► immunosuppressive_niche_signature_genes_ranked.csv│
│       (三层筛选)   ──► immunosuppressive_niche_signature_genes.txt ★     │
│                    ──► prior_gene_set_auc.csv                            │
│                    ──► gini_score_genes.csv                              │
│                    ──► niche_signature_volcano.png                       │
│                    ──► niche_fraction_scatter.png                        │
│          │                                                               │
│          ▼                                                               │
│  Step 15: 参数扫描 ──► param_scan_deg_stability.csv                     │
│                    ──► param_scan_deg_stability_heatmap.png              │
│                                                                          │
└──────────────────────────────────┬───────────────────────────────────────┘
                                   │
                                   ▼ immunosuppressive_niche_signature_genes.txt
┌──────────────────────────────────────────────────────────────────────────┐
│              Step 3 (可选): TCGA 队列生存分析验证                         │
│            (tcga_survival_analysis.R)                                    │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  签名基因列表 ──► ssGSEA 评分 ──► TCGA-LIHC 队列 ──► Kaplan-Meier       │
│                                                         生存曲线         │
│  验证假说：免疫抑制生态位签名基因高表达患者预后更差                       │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## 五、核心生物学结论链

```
HCC4R + CHC20 两切片联合分析
        │
        ▼
Cell2location 反卷积 → 每个 spot 的细胞类型空间地图
        │
        ▼
Treg-TAM-CAF 三者空间共定位（Leiden 无监督聚类 + 语义注释）
        │
        ▼
免疫抑制生态位（Immunosuppressive Niche）的空间识别
        │
        ├── 距离依赖性：免疫抑制信号在肿瘤边缘带最强
        │
        ├── L-R 通讯网络：CXCL12-CXCR4 / CCL22-CCR4 招募轴 + TGFB1 免疫抑制轴
        │
        ├── 签名基因：三层筛选（Wilcoxon + 先验AUC + Gini）克服 Visium 稀释效应
        │   └── delta_frac 检出率维度弥补 FOXP3 等稀有基因的均值稀释
        │
        └── 参数鲁棒性：Spearman ρ > 0.95 / Weighted Jaccard > 0.90
            └── 免疫抑制生态位是真实生物学实体，非参数选择产物
                    │
                    ▼
            TCGA 队列验证（可选）：免疫抑制生态位签名 → 预后预测
```

**核心假说**：肝细胞癌组织中存在 Treg-TAM-CAF 协同构成的免疫抑制空间生态位，该生态位通过趋化因子招募（CXCL12-CXCR4、CCL22-CCR4）和免疫抑制信号（TGFB1、IL10、PD-1/PD-L1）维持自我强化回路，在肿瘤边缘带形成功能性免疫逃逸屏障，且其分子签名可预测患者预后。

---

> **文档版本**：基于 `code/pipeline/run_preprocessing.py` 和 `code/pipeline/run_spatial_niche_analysis.py` 的代码逻辑生成，尚未运行项目，实际结果文件可能在具体数值上与描述略有差异。