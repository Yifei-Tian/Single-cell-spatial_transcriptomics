# 分析流程输出文件说明

本文档对本项目所有 `results/` 目录下生成的输出文件进行逐一解释，说明每个文件的内容格式、如何阅读和解读，以及其背后揭示的生物学意义。

**当前分析模式：HCC4R 主分析 + CHC20 独立验证**
- **HCC4R**（Discovery Cohort）：主分析切片，运行 Step 1 + Step 2 + Step 3，结果保存在 `results/HCC4R/`
- **CHC20**（Validation Cohort）：独立验证切片，独立运行 Step 1 + Step 2，结果保存在 `results/CHC20/`，用于 Figure 5 跨队列验证
- 运行入口：`python code/run_hcc4r.py`（Step 1 + Step 2 + 可选 Step 3）
- 论文图表：`python code/run_hcc4r.py --paper-figures` 或 `--figures-only`

分析流程共三个核心步骤，输出文件按步骤组织。

---

## 一、Step 1 — 数据预处理与 Cell2location 反卷积（`pipeline/run_preprocessing.py`）

> 输出根目录：`results/HCC4R/`（主分析），`results/CHC20/`（验证分析结构相同）

### 1.1 `adata_vis_post.h5ad`

**文件格式**：HDF5 格式的 AnnData 对象（单细胞/空间组学标准存储格式）

**内容说明**：

这是 HCC4R 切片经过 Cell2location 贝叶斯反卷积后的核心数据文件，是后续所有空间分析的起点。文件内部包含：

- `adata.X`：每个 spot 的原始 count 表达矩阵（行 = spot，列 = 基因）
- `adata.obs`：每个 spot 的元数据（spot barcode、组织标记等）
- `adata.obsm["spatial"]`：每个 spot 的二维空间坐标（来自 Visium 芯片物理布局）
- `adata.obsm["means_cell_abundance_w_sf"]`：Cell2location 后验估计的细胞类型丰度矩阵（行 = spot，列 = 细胞类型），这是反卷积的主要输出

**如何理解**：

10x Visium 空间转录组每个 spot 直径约 55 μm，每个 spot 捕获的 RNA 来自多种细胞类型的混合信号。Cell2location 是一种贝叶斯层次模型，利用 scRNA-seq 参考数据中各细胞类型的基因表达特征（参考签名），通过负二项分布似然函数将每个 spot 的混合信号"反卷积"为各细胞类型的估计数量（后验均值 `means_cell_abundance_w_sf`）。

**揭示的生物学现象**：

反卷积结果直接回答了"每个空间位置存在哪些细胞类型"这一核心问题。通过观察 `means_cell_abundance_w_sf` 矩阵，可以看到：
- 肿瘤实质（Hepatocyte 丰度高）在组织切片中的分布
- 免疫细胞（Treg、T/NK、Myeloid）的浸润区域及其空间梯度
- 基质细胞（Fibroblast）形成的物理屏障位置

这是从"基因表达"到"细胞类型空间分布"的关键转化步骤。

---

### 1.2 `adata_sc_post.h5ad`

**文件格式**：HDF5 格式的 AnnData 对象

**内容说明**：

scRNA-seq 参考数据集经过 Cell2location RegressionModel 训练后的 AnnData，包含：

- 原始 scRNA-seq 表达矩阵（基础质控后）
- `adata.uns["mod"]`：RegressionModel 训练后的模型参数（后验推断结果）
- 细胞类型标签（obs 列）和 Treg 重注释结果

**如何理解**：

这是 Cell2location 的中间产物，存储了从 scRNA-seq 学习到的每种细胞类型的"参考签名矩阵"。该文件可用于：① 检查哪些基因对各细胞类型区分贡献最大；② 为其他下游分析（如验证签名质量）提供参考。

---

### 1.3 `spot_cell_proportion_HCC4R.csv`

**文件格式**：CSV 表格

**列说明**（每列代表一种细胞类型，每行代表一个 spot）：

| 列名 | 含义 |
|------|------|
| `spot_id` | Visium spot 的唯一标识符（如 `ACGCCTGACACGCGCT-1`） |
| `Hepatocyte` | 该 spot 中肝细胞比例（0–1） |
| `Malignant` | 恶性肿瘤细胞比例 |
| `Treg` | 调节性 T 细胞比例 |
| `T/NK` | T 细胞和 NK 细胞的合并比例 |
| `Myeloid` | 髓系细胞（包括 TAM、MDSC 等）比例 |
| `Fibroblast` | 成纤维细胞（CAF）比例 |
| `HSC` | 肝星状细胞比例 |
| `B cell` | B 细胞比例 |
| `Endothelial` | 内皮细胞比例 |
| ... | 其他细胞类型 |

所有细胞类型比例之和为 1（归一化后）。

**如何理解**：

这是将 `adata_vis_post.h5ad` 中的绝对丰度值归一化为相对比例的结果表。每行是一个 spot 的细胞组成"快照"，直接反映该空间位置的局部微环境细胞构成。

**揭示的生物学现象**：

- `Treg` 列的数值分布可揭示肿瘤组织中免疫抑制压力的空间分布格局
- `Hepatocyte` 高值区域对应肿瘤实质（肿瘤核心区），低值区域对应间质/免疫浸润区
- `Treg` 与 `Myeloid` 同时高值的 spot 提示可能的免疫抑制微生态位区域

---

### 1.4 `shared_genes_HCC4R.txt`

**文件格式**：纯文本，每行一个基因名

**内容说明**：scRNA-seq 参考数据与 HCC4R Visium 切片**共同基因列表**，即最终用于 Cell2location 建模的基因集。

**如何理解**：交集（scRNA ∩ HCC4R）保证了反卷积中参考签名与空间数据使用相同的基因空间。共有基因数量通常在 2,000–5,000 个之间（取决于测序深度和 Visium 捕获效率）。

---

### 1.5 `scrna_tsne_celltype.png`（论文 Figure 1B）

**文件格式**：PNG 图像

**内容说明**：

scRNA-seq 单细胞数据的 tSNE（或 UMAP）嵌入可视化图，每个点代表一个细胞，颜色按细胞类型分配（高饱和度离散色板），并在每种细胞类型聚集的中心位置标注类型名称，右侧图例完整展示颜色-类型对应关系，坐标轴以箭头样式标注 tSNE 1 / tSNE 2 方向。

**如何理解**：

此图是分析流程的起点可视化，展示了 scRNA-seq 参考数据集中各细胞类型在低维嵌入空间中的分布格局。不同细胞类型的点云应形成清晰分离的簇，且生物学相近的细胞类型（如 Treg 与 T/NK）在空间上相邻。与参考论文 `docs/plots/1.png` 风格一致。

**揭示的生物学现象**：

- 各细胞类型的分离程度反映参考数据集的质量，分离良好的聚类说明 scRNA-seq 数据捕获了足够的转录异质性
- Treg 细胞应与其他 T/NK 细胞聚簇相邻但明显分开（若分辨率足够），这验证了 Treg 亚群鉴定的准确性
- Myeloid、Fibroblast 等基质细胞与免疫细胞（T/NK、Treg）的分离反映了肝癌组织微环境中细胞类型组成的多样性

---

### 1.6 `scrna_celltype_marker_heatmap.png`（论文 Figure 1C）

**文件格式**：PNG 图像

**内容说明**：

细胞类型 Marker 基因平均表达热图，X 轴为细胞类型，Y 轴为各细胞类型的特异性标志基因，颜色采用 Yellow-Black-Purple 三色渐变（paper_ybp 配色：亮黄=高表达，黑色=中等表达，深紫=低/无表达），值为 CP10K 归一化 + log1p 变换后各细胞类型内的均值，再经跨细胞类型 Z-score 标准化（便于横向对比）。与参考论文 `docs/plots/2.png` 风格一致。

**主要 Marker 基因组**：

| 细胞类型 | 关键 Marker 基因 |
|---------|----------------|
| T cell | IL7R, CD3G, CD2, ITM2A, CD3D |
| Myeloid | LYZ, AIF1, RNASE1, C1QB, HLA-DRA |
| NK | GNLY, GZMB, KLRD1, KLRF1 |
| B cell | B3GNT7, MS4A1, BANK1, CD79A |
| Malignant | APOA2, ALB, APOA1, AMBP, APOH |
| Endothelial | PECAM1, CDH5, SPARCL1, STC1 |
| Epithelial | INSR, KRT18, KRT19, EPCAM |
| Plasma cell | JCHAIN, MZB1, IGLL5 |
| HSC | RGS5, COL1A1, ACTA2, PDGFRB |

**如何理解**：

热图中每列代表一种细胞类型，理想情况下各类型的特异性 Marker 在对应列中呈现亮黄色（高 Z-score），而在其他列中保持深紫色（低表达）。这种"对角线"高亮模式验证了细胞类型注释的可靠性。

**揭示的生物学现象**：

- ALB、APOA2 等白蛋白相关基因在 Malignant/Hepatocyte 类中特异性高表达，确认了肿瘤细胞的肝细胞来源
- FOXP3 在 Treg 簇的检出模式（即使表达量绝对值低，Z-score 相对差异仍清晰可见）验证了调节性 T 细胞的鉴定
- COL1A1、ACTA2 在 HSC/Fibroblast 中的高表达揭示了肝星状细胞激活后的基质重塑能力，与免疫抑制微环境的 CAF 功能相关

---

### 1.7 `t_cell_dotplot_horizontal.png`（论文 Figure 1D）

**文件格式**：PNG 图像

**内容说明**：

横版气泡图（Dot Plot），展示 T/NK 细胞各子簇（subclusters，分辨率 res.3 下的 Leiden 簇）中四个关键标志基因的表达情况：

**CD3D（泛 T 细胞标志物）：**

- **作用：** 这是一个门槛。所有真正的 T 细胞都应该表达它。
- **看图：** 你可以看到这 4 个基因中，CD3D 在绝大多数簇（1, 11, 13 等）都有较大且偏蓝的气泡，证明这些基本上都是 T 细胞。

**CD4（辅助性 T 细胞和treg均表达）：**

- **作用：** Treg 属于 CD4+ T 细胞的一个亚群。所以 Treg 细胞所在的簇，CD4 也必须有一定程度的表达。

**FOXP3（Treg 的核心转录因子，特异性最强）：**

- **作用：** 这是 Treg 细胞最核心的主调控转录因子。在 T 细胞中，**只有 Treg 才会特异性地表达 FOXP3**。这是你破局的关键点。

**IL2RA / CD25（Treg 活化标志，在稳定 Treg 中高表达）：**

- **作用：** 稳定的 Treg 细胞表面会高表达 IL2RA。它通常与 FOXP3 协同出现。

气泡大小表示在某一个特定的细胞簇（x轴的数字）中，**有多少比例的细胞**表达了该基因。气泡越大，说明该群体中表达这个基因的细胞占比越高。

颜色深浅表示在那些表达了该基因的细胞中，**基因表达量的平均水平有多高**。颜色越蓝、越深，说明该基因的表达量极高；颜色越白、越浅，说明表达量较低。

**揭示的生物学现象**：

FOXP3 是 Treg 的"主调控转录因子"（master transcription factor），其特异性表达是鉴定 Treg 的金标准。图中清晰显示 FOXP3 仅在某一特定簇（如簇9）中高表达，表明该簇为真正的调节性 T 细胞亚群，而非普通辅助 T 细胞或效应 T 细胞，为后续空间 Treg 分析奠定细胞类型定义基础。

---

### 1.8 `regression_training_history_HCC4R.png`

**文件格式**：PNG 图像

**内容说明**：

Cell2location 第一阶段 RegressionModel（单细胞参考签名学习）的训练损失曲线，x 轴为训练 epoch 数，y 轴为 ELBO（证据下界，Evidence Lower BOund）。

**如何理解**：

ELBO 是变分推断（Variational Inference）中的优化目标，其绝对值越大（负 ELBO 越小）说明模型对数据的拟合越好。曲线应在若干 epoch 后趋于平稳（收敛），若曲线仍在下降则说明训练轮次不足。

**揭示的生物学现象**：

收敛的训练曲线保证了细胞类型参考签名的可靠性——只有当 RegressionModel 充分学习到每种细胞类型的特异性基因表达模式后，下一阶段的空间反卷积才能准确区分不同细胞类型。

---

### 1.9 `cross_slice_comparison/` 目录（可选，多切片一致性验证图）

当 Step 1 同时传入 HCC1R 配对验证切片（`--no-sample2` 参数未启用时）生成。

#### `cross_slice_mean_proportion_comparison.png`

**内容**：HCC4R 与 HCC1R 各细胞类型**全切片平均比例**的并排条形图。

**如何理解**：若两个切片的细胞类型组成比例趋势相近（如 Hepatocyte 均为最高占比，Treg 均偏低但存在），则说明两个配对样本的整体微环境构成具有代表性，支持一致性假设。

#### `cross_slice_treg_distribution.png`

**内容**：HCC4R 与 HCC1R Treg 比例分布的核密度估计（KDE）对比曲线。

**揭示的生物学现象**：Treg 在肝癌组织中的浸润通常是稀疏但局灶性富集的（分布呈右尾），这种"局灶性免疫抑制"模式的跨样本一致性是免疫逃逸机制普遍存在的空间证据。

#### `cross_slice_celltype_boxplot.png`

**内容**：Treg、Myeloid、Fibroblast、Hepatocyte 四种关键细胞类型在两切片中比例分布的箱线图对比。

---

## 二、Step 2 — 空间免疫抑制生态位分析（`pipeline/run_spatial_niche_analysis.py`）

> 输出根目录：`results/HCC4R/spatial_niche/`（主分析）

**关键参数**：kNN 邻居数 k=15，niche_high 分位数阈值 q=0.85（均经参数扫描验证为最优值，详见 `param_scan_deg_stability_heatmap.png`）

### 2.1 `spatial_niche/spatial_niche_scores.csv`

**文件格式**：CSV 表格（行 = spot，列 = 多种评分和标签）

**关键列说明**：

| 列名 | 类型 | 含义 |
|------|------|------|
| `spot_id` | 字符串 | Spot 唯一标识符（Visium barcode） |
| `spatial_x` / `spatial_y` | 浮点数 | 空间坐标（像素单位） |
| `Hepatocyte` / `Malignant` / `Treg` / `Myeloid` / `Fibroblast` / `T/NK` / `HSC` 等 | 浮点数（0–1） | 各细胞类型归一化比例 |
| `immunosuppressive_gene_score` | 浮点数 | 免疫抑制标志基因（FOXP3、IL2RA、CTLA4、TIGIT、LAG3、TGFB1、IL10 等）的 CP10K+log1p 归一化平均表达量 |
| `neighborhood_cluster` | 字符串（如"0","1",...） | Leiden 邻域组成聚类分配的 cluster 编号（k=15 kNN，resolution=0.5） |
| `niche_semantic_label` | 字符串 | 语义注释标签：`immunosuppressive_niche` 或 `cluster_X` |
| `Treg_like_score` | 浮点数 | Z-score(Treg比例) + Z-score(免疫抑制基因评分) 之和 |
| `immune_stroma_score` | 浮点数 | Treg + T/NK + Myeloid + Fibroblast 四项 Z-score 之和 |
| `immunosuppressive_niche_score` | 浮点数 | 综合免疫抑制 niche 评分（5项 Z-score 加权求和，含邻域 Hepatocyte 均值） |
| `niche_high` | 布尔值 | 该 spot 是否属于高免疫抑制区域（评分 ≥ **85 百分位**，经参数扫描验证的最优阈值） |
| `spatial_region` | 字符串 | 规则化空间区域标注（`tumor_core`/`tumor_edge`/`stroma_immune`/`other`） |
| `hep_high` | 布尔值 | 该 spot 是否属于高 Hepatocyte 比例区域（超过 75 百分位） |
| `distance_to_hep_high` | 浮点数 | 到最近高肝细胞区域的欧式距离 |
| `neighbor_Hepatocyte` 等 | 浮点数 | 空间半径邻域内各细胞类型比例的均值（平滑后的局部微环境信息） |

**如何理解整个评分体系**：

本分析采用两阶段策略识别免疫抑制生态位：

**阶段 1：数据驱动的无监督发现（划定"社区"边界）**

- 动作 1：获取邻域信息（平滑化）。孤立地看一个 spot 是不准的。这一步计算了每个 spot 周围最近的 15 个"邻居"（**k=15，参数扫描验证的最优值**）的细胞组成均值，消除了单点噪声，真正抓取了"微环境"特征。
- 动作 2：Leiden 聚类。算法根据上面算出的"街道平均特征"进行无监督聚类（resolution=0.5），把细胞组成相似的区域归为一类，生成了 `neighborhood_cluster`（例如 Cluster 0, 1, 2...）。

**阶段 2：功能评分与语义注释（给"社区"贴标签）**

- 动作 3：多维特征考核。针对每一个聚类（Cluster），计算四个关键指标的均值：**Treg 比例、Myeloid 比例（通常代表肿瘤相关巨噬细胞 TAM）、Fibroblast 比例（肿瘤相关成纤维细胞 CAF）、免疫抑制基因评分**。
- 动作 4：最高分胜出（命名为 immunosuppressive_niche）。在这四个维度上综合得分最高的那个 Cluster，被正式赋予了 `niche_semantic_label = "immunosuppressive_niche"` 的头衔。

**附加阶段：连续评分与高低分组（用于下游差异分析）**

- 动作 5：综合打分（Score）。为**每一个** spot 计算连续的 `immunosuppressive_niche_score`（5项 Z-score 的加权和）。
- 动作 6：切分阈值（High vs Low）。通过设定阈值（**超过 85 百分位，经参数扫描验证的最优阈值**），硬性划分出 `niche_high` 区域，用于签名基因提取的分组。

**揭示的生物学现象**：

- `niche_semantic_label = "immunosuppressive_niche"` 的 spot 集合定义了**肿瘤免疫抑制微生态位**的空间范围，即 Treg、TAM（Myeloid）和 CAF（Fibroblast）三类细胞同时局部富集的区域
- 这类区域在生物学上对应肿瘤逃逸的"免疫沙漠"或"免疫排除"区域
- `spatial_region` 的分层（肿瘤核心→肿瘤边缘→免疫基质区域）反映了肿瘤组织内部微环境的空间异质性

---

### 2.2 `spatial_niche/spatial_niche_parameters.csv`

**文件格式**：CSV 表格（两列：`parameter` 和 `value`）

**内容说明**：记录本次分析使用的所有关键参数及运行时计算得到的统计量，包括：

| 参数 | 含义 |
|------|------|
| `n_neighbors_knn` | 邻域组成聚类时的 k-NN 邻居数（**当前值：15**，经参数扫描验证） |
| `leiden_resolution` | Leiden 聚类分辨率（默认0.5） |
| `n_neighborhood_clusters` | 实际识别出的 Leiden 聚类数量 |
| `neighbor_radius_multiplier` | 空间半径邻域倍增系数（默认1.25） |
| `neighbor_radius` | 实际使用的搜索半径（单位：坐标像素） |
| `hep_high_threshold` | 高肝细胞区域的 Hepatocyte 比例阈值（75 百分位值） |
| `niche_high_threshold` | niche_high 评分阈值（**85 百分位值，经参数扫描验证**） |
| `n_niche_high` | 被标记为 niche_high 的 spot 总数 |
| `niche_high_quantile` | niche_high 分位数设置（**当前值：0.85**） |
| `marker_genes_used` | 实际用于计算免疫抑制基因评分的基因列表 |
| `per_sample_neighbors` | 是否启用按切片内建立空间邻域（联合分析时为 True，单切片分析时为 False） |

**如何理解**：

该文件是分析的"元数据记录"，作用包括：① 重现性保证（记录参数可供他人完整复现）；② 质控检查（如 `n_neighborhood_clusters` 过多或过少均提示聚类参数需调整）；③ 报告依据（论文方法部分的参数来源）。

---

### 2.3 `spatial_niche/sensitivity_analysis.csv`

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| `param_mode` | 参数类型（`knn` 或 `radius`） |
| `param_label` | 参数具体设置（如 `kNN k=10`、`radius×1.0`） |
| `spearman_rho` | Spearman 秩相关系数（与参考配置 k=15 相比的评分排序一致性） |
| `weighted_jaccard` | 加权 Jaccard 相似度（连续版本，对高分区域天然加权） |
| `is_reference` | 是否为参考配置（k=15 kNN）|

**如何理解**：

敏感性分析（v2 连续指标版本）扫描了 kNN 邻居数（k=10/15/20）和空间半径倍增系数（×1.0/1.25/1.5）共 6 种参数组合，采用两个连续性稳健指标衡量结果稳定性：
- `spearman_rho`：衡量全局评分排序的一致性，不受硬截断影响；
- `weighted_jaccard`：基于 Min-Max 归一化连续向量计算，对高分区域天然加权。

关键判断标准：若 `spearman_rho > 0.90` 且 `weighted_jaccard > 0.80`，则说明 niche 识别结果**对参数选择不敏感**，结论具有稳健性。

**揭示的生物学现象**：

稳健的 niche 边界意味着免疫抑制微生态位是一种真实存在的、结构性的空间特征，而非因任意选定某个阈值造成的人为划分。这增强了该生态位作为真实生物学实体的可信度。

---

### 2.4 `spatial_niche/neighborhood_cluster_stats.csv`

**文件格式**：CSV 表格（行 = Leiden cluster，列 = 各维度统计量）

**列说明**：

| 列名 | 含义 |
|------|------|
| `mean_treg` | 该 cluster 内所有 spot 的 Treg 比例均值 |
| `mean_myeloid` | Myeloid 比例均值 |
| `mean_fibroblast` | Fibroblast 比例均值 |
| `mean_immune_score` | 免疫抑制基因评分均值 |
| `rank_mean_treg` 等 | 各维度的降序排名（排名越小说明该 cluster 该维度越高） |
| `combined_rank` | 四维排名之和（越小说明综合免疫抑制特征越强） |

**如何理解**：

`combined_rank` 最小的 cluster 被注释为 `immunosuppressive_niche`，这是一种基于多维证据权衡的无偏注释方法，避免了仅依赖单一指标（如只看 Treg）导致的偏差。

**揭示的生物学现象**：

该表量化展示了"免疫抑制生态位 cluster"与其他 cluster 在免疫抑制特征上的差异程度。若目标 cluster 的 `combined_rank` 远低于次低 cluster，则说明免疫抑制信号的空间集中性非常突出，具有明显的"hot spot"特征。

---

### 2.5 `spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv`

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| `gene` | 基因名（HGNC 符号） |
| `mean_high` | 该基因在 `niche_high`（≥85百分位）spot 中的平均 log1p 表达量 |
| `mean_low` | 该基因在 `niche_low` spot 中的平均 log1p 表达量 |
| `log2_fc` | log₂(mean_high+1) - log₂(mean_low+1)，衡量在免疫抑制区域的富集倍数 |
| `frac_high` | `niche_high` 组中检出该基因（表达>0）的 spot 比例 |
| `frac_low` | `niche_low` 组中检出该基因（表达>0）的 spot 比例 |
| `delta_frac` | `frac_high - frac_low`：检出率差值（弥补均值 log2FC 对稀有基因不敏感的缺陷） |
| `composite_score` | 综合排序分数（0.5×norm(log2FC) + 0.5×norm(delta_frac)） |
| `pvalue` | Wilcoxon 秩和检验 p 值 |
| `fdr` | BH-FDR 校正后的 p 值 |

按 `composite_score` 降序排列。

**如何理解**：

本流程采用**三层筛选策略**：
- **Layer 1**（Wilcoxon + FDR + 综合排序）：筛选条件 FDR<0.05 且（log2FC>0.5 OR delta_frac>0.10）
- **Layer 2**（先验功能基因集 AUC 检验）：Treg/TAM/CAF 三组先验基因中 AUC>0.6 且 delta_frac>0.05 的强制合并
- **Layer 3**（Gini Index 特异性评分）：Gini>0.3 且 log2FC>0，用于检测 FOXP3 等局灶性高表达稀有免疫基因

`delta_frac` 是针对 Visium spot-level 数据的关键补充维度——由于每个 spot 覆盖多个细胞，FOXP3 等稀有 Treg 标志基因因细胞稀释效应在大多数 spot 中为 0，仅用均值 log2FC 无法检出；`delta_frac` 度量"niche_high 组有更多 spot 能检出该基因"，弥补了这一缺陷。

**揭示的生物学现象**：

签名基因集代表了肿瘤局部免疫抑制环境的分子特征。典型上调基因通常包括：
- **Treg 标志基因**：FOXP3、IL2RA、CTLA4、TIGIT（免疫检查点分子）
- **免疫抑制细胞因子**：TGFB1、IL10（直接抑制效应 T 细胞活性）
- **趋化因子**：CXCL12、CCL22（招募 Treg 入侵肿瘤区域的信号分子）
- **基质重塑相关基因**：FAP、ACTA2（CAF 激活标志，参与物理屏障构建）

---

### 2.6 `spatial_niche/immunosuppressive_niche_signature_genes.txt`

**文件格式**：纯文本，每行一个基因名

**内容**：`immunosuppressive_niche_signature_genes_ranked.csv` 中基因列的纯文本版本，是传递给 `tcga_survival_analysis.R` 进行 TCGA 队列验证的**关键接口文件**。

---

### 2.7 `spatial_niche/prior_gene_set_auc.csv`

**文件格式**：CSV 表格

**内容**：Layer 2 先验基因集（Treg/TAM/CAF 三组）的 AUC 检验结果，包含 `frac_high`、`frac_low`、`delta_frac` 列。

**如何理解**：AUC > 0.6 且 delta_frac > 0.05 的先验基因被强制纳入签名基因列表，即使其 Wilcoxon p 值未通过 FDR 截断（因稀释效应导致统计功效不足），保证了生物学上已知重要的免疫调控基因不被漏掉。

---

### 2.8 `spatial_niche/gini_score_genes.csv`

**文件格式**：CSV 表格

**内容**：Layer 3 Gini Index 特异性评分筛选出的基因列表（Gini > 0.3 且 log2FC > 0）。

**如何理解**：Gini 指数衡量基因表达在 spot 间的"局灶性"——Gini 高意味着该基因仅在少数 spot 中高表达（局灶富集），这正是 FOXP3 等稀有免疫细胞标志基因的典型模式。传统的均值比较无法检出此类基因，Gini 指数是专门针对空间局灶性信号设计的补充指标。阈值从原来的 0.5 降低到 0.3，覆盖更多稀有免疫基因。

---

### 2.9 `spatial_niche/param_scan_deg_stability.csv`

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| `k` | kNN 邻居数（8 / 10 / 15 / 20 / 25） |
| `quantile` | niche_high 分位数阈值（0.70 / 0.75 / 0.80 / 0.85） |
| `n_sig_deg` | 该参数组合下显著 DEG 数量（FDR<0.05 且 log2FC>0.5） |
| `mean_log2fc_topN` | Top-N 基因的平均 log2FC |
| `top_genes_str` | Top 基因名列表（字符串形式） |

**如何理解**：

k × quantile 二维网格搜索结果原始数据，是绘制 `param_scan_deg_stability_heatmap.png` 的数据来源。**当前最优参数 k=15, q=0.85 即来自此表的分析结果**——在该参数组合处，显著 DEG 数量处于"高原区"（与相邻格子结果相近），表明参数选择稳健。

---

### 2.10 `spatial_niche/plots/` 可视化图表

#### `spatial_hepatocyte.png`

**内容**：Hepatocyte 细胞比例在切片空间坐标中的分布，使用 Yellow-Black-Purple 渐变色标（paper_ybp 配色：深紫=低，黑色=中，亮黄=高）。

**如何理解**：高 Hepatocyte 比例区域（亮黄色）对应肿瘤实质，即肝癌细胞密集区；低比例区域（深紫色）对应肿瘤间质（免疫浸润区、基质区）。此图是定义肿瘤核心区域的空间依据。

#### `spatial_treg.png`

**内容**：Treg 细胞比例的空间分布图（paper_ybp 色标，亮黄色区域为 Treg 富集热点）。

**如何理解**：Treg 在肝癌中的分布通常并非均匀弥散，而是局灶性富集在特定区域（如肿瘤边界或免疫浸润区）。识别 Treg 富集的"热点"区域是本研究的核心目标之一。

**揭示的生物学现象**：局灶性 Treg 富集区域通常对应肿瘤内部的免疫抑制微生态位，这些区域中效应 T 细胞的杀伤功能受到最强烈的抑制，是肿瘤免疫逃逸的关键空间节点。

#### `spatial_myeloid.png` / `spatial_fibroblast.png`

**内容**：Myeloid（髓系细胞）和 Fibroblast（成纤维细胞）比例的空间分布（paper_ybp 色标）。

**揭示的生物学现象**：
- Myeloid 细胞（包括肿瘤相关巨噬细胞 TAM）通常与 Treg 共分布，构成 Treg 招募与维持的细胞环境（CCL22-CCR4 轴）
- Fibroblast 高比例区域对应肿瘤基质，CAF 通过分泌 CXCL12 等趋化因子和重塑细胞外基质参与免疫抑制

#### `spatial_immunosuppressive_niche_score.png`

**内容**：综合免疫抑制 niche 评分的空间分布，使用 Yellow-Black-Purple 渐变色标（paper_ybp 配色：深紫=低分，黑色=中等，亮黄=高分热点）。

**如何理解**：paper_ybp 色标的选择使高评分区域（免疫抑制热点）在视觉上非常突出（亮黄色）。整个切片的"评分地图"直观揭示了免疫抑制压力的空间分布格局。

**揭示的生物学现象**：高评分区域即免疫抑制微生态位的空间轮廓，该区域通常位于肿瘤实质与免疫浸润区的交界地带（即免疫细胞进入肿瘤的"前哨区域"），是制定靶向免疫治疗策略时需要重点关注的空间区域。

#### `spatial_neighborhood_clusters.png`

**内容**：Leiden 邻域组成聚类结果的空间分布图，每个 cluster 用不同颜色标识（高饱和度离散色板，仿论文风格）。

**如何理解**：这是分析方法的"中间过程可视化"，展示无监督聚类对肿瘤组织空间结构的划分。空间上相邻且颜色相同的 spot 组成一个具有相似局部细胞组成的"邻域"。

**揭示的生物学现象**：聚类结果反映了组织的**空间层级结构**：不同 cluster 对应肿瘤的不同功能区域（如肿瘤核心、免疫浸润区、基质区等）。这种无监督发现的空间分区与病理学上的组织学分区往往高度一致。

#### `spatial_niche_semantic_labels.png`

**内容**：语义注释后的 niche 标签空间分布图，其中紫色（`#9467bd`）标注 `immunosuppressive_niche`，灰色标注其他 cluster。

**如何理解**：这是分析的核心结论图之一。紫色区域直接标识了经数据驱动方法识别的免疫抑制微生态位的空间轮廓和面积占比。

#### `spatial_region_labels.png`

**内容**：规则化空间区域标注图，用四种颜色标注：`tumor_core`（红色）、`tumor_edge`（橙色）、`stroma_immune`（蓝色）、`other`（灰色）。

**如何理解**：

- **tumor_core**：Hepatocyte 比例 > 75 百分位，即肿瘤实质区
- **tumor_edge**：非高肝细胞区域但与 tumor_core 相邻，且免疫基质评分较高，即肿瘤侵袭边界
- **stroma_immune**：免疫基质评分 > 60 百分位的非肿瘤实质区，即富含免疫细胞和基质细胞的区域

**揭示的生物学现象**：肿瘤的"边界区域"（tumor_edge）通常是免疫细胞浸润最活跃但同时免疫抑制最强烈的地带，是肿瘤-免疫细胞相互作用的主要战场。

#### `spatial_niche_high_score_spots.png`

**内容**：niche_high 二值分布图，标记 immunosuppressive_niche_score ≥ **85 百分位**（q=0.85，经参数扫描验证的最优阈值）的 spot 空间位置。

**如何理解**：相比语义标签图（基于 Leiden 聚类），该图是基于连续评分硬截断的结果，两者并排比较（见下一图）可验证两种方法的一致性。

#### `spatial_niche_cluster_vs_score_comparison.png`

**内容**：聚类法（`niche_semantic_label`，左图）与评分法（`niche_high` q=0.85，右图）并排对比图，同一切片两种方法的结果可视化。

**如何理解**：两种方法的空间分布高度重叠则验证了两种独立发现方法的一致性，进一步确认免疫抑制生态位的真实性。若仅一种方法发现的区域，可能代表两种方法各自捕获到的不同侧面。

#### `hepatocyte_vs_treg_niche_score.png`

**内容**：Hepatocyte 比例（x 轴）vs Treg 比例（y 轴）散点图，颜色编码综合 niche 评分（magma 色标），虚线标注 75 百分位分界线划分四象限。

**如何理解**：重点关注右上象限（Hepatocyte 高 + Treg 高）中颜色最深的点，这些 spot 是肿瘤实质与 Treg 共定位最明显的区域。

**揭示的生物学现象**：免疫细胞（包括 Treg）无法突破物理屏障，只能在肿瘤核心（高 Hepatocyte 区）的外围大量聚集。这通常与肿瘤外围致密的细胞外基质有关。

#### `distance_to_hep_high_vs_niche_score.png`

**内容**：距离高肝细胞区域（肿瘤核心）的距离（x 轴，等频分箱）vs 平均免疫抑制 niche 评分（y 轴）的折线图，附带标准误差棒。

**如何理解**：若折线呈现"距肿瘤核心越近，免疫抑制评分越高"的趋势，则说明免疫抑制信号具有**距离依赖性**（肿瘤-近端区域受到更强的免疫抑制压力）。

**揭示的生物学现象**：肿瘤实质边界区域的强免疫抑制信号与肿瘤招募 Treg 和分泌趋化因子（如 CXCL12）的机制一致，是肿瘤主动塑造免疫抑制微环境的空间证据。

#### `region_score_boxplots.png`

**内容**：四个空间区域（`tumor_core`、`tumor_edge`、`stroma_immune`、`other`）中 `Treg_like_score`（左图）和 `immunosuppressive_niche_score`（右图）的箱线图。

**如何理解**：若 `stroma_immune` 或 `tumor_edge` 区域的箱体中位线最高，说明免疫抑制核心集中在免疫浸润区而非纯肿瘤实质，提示 Treg 的主要作用位点在"肿瘤-免疫交界面"。

#### `celltype_niche_correlation.png`

**内容**：各细胞类型比例与 niche 评分之间的 Pearson 相关系数热图（vlag 配色：红色=正相关，蓝色=负相关）。

**如何理解**：关注 Treg 与 Myeloid 之间的相关系数、以及这两者与 Hepatocyte 之间的相关系数。

**揭示的生物学现象**：
- **Treg 与 Myeloid 正相关**：支持 TAM 通过分泌 CCL22 招募 Treg 的机制假设（TAM→Treg 信号轴）
- **Treg 与 Hepatocyte 正相关或负相关**：揭示肿瘤细胞与 Treg 的空间关联模式
- **Fibroblast 与 Myeloid 正相关**：提示 CAF-TAM 共定位在基质区形成双重免疫屏障

#### `lr_communication_heatmap.png`

**内容**：配体-受体（L-R）通讯分析热图，共两个子图：

- **左图**：关键 L-R 信号对（行）在 `Niche-High` vs `Niche-Low` 两组（列）中的归一化信号强度热图（YlOrRd 色标）
- **右图**：每个 L-R 对的 log₂FC（Niche-High/Niche-Low），红色=在免疫抑制 niche 中上调，蓝色=下调

**关键 L-R 信号对及其生物学意义**：

| L-R 对 | 信号通路 | 生物学功能 |
|--------|----------|-----------|
| CXCL12–CXCR4 | SDF-1 轴 | TAM 和基质细胞分泌 CXCL12，招募表达 CXCR4 的 Treg 进入肿瘤，是肝癌 Treg 招募的核心机制 |
| CCL22–CCR4 | 趋化因子轴 | M2 巨噬细胞分泌 CCL22，结合 Treg 表面 CCR4 受体，实现定向招募 |
| TGFB1–TGFBR1 | TGF-β 轴 | TGF-β 是 Treg 产生和维持的关键细胞因子，同时直接抑制效应 T 细胞 |
| PD-1–PD-L1 | 免疫检查点 | 肿瘤细胞/Myeloid 细胞上的 PD-L1 结合 T 细胞表面 PD-1，直接触发效应 T 细胞功能失活 |
| IL10–IL10RA | 抗炎细胞因子轴 | Treg 和 TAM 分泌 IL-10，通过 IL10RA 信号在局部造成免疫耐受环境 |
| TIGIT–NECTIN2 | 免疫检查点 | Treg 和耗竭 T 细胞表达 TIGIT，与 APC 上的 NECTIN2 结合触发抑制信号 |
| LAG3–MHC-II | 免疫检查点 | LAG3 结合 MHC-II 分子，与 TIGIT/PD-1 协同参与 T 细胞耗竭 |
| SPP1–CD44 | 骨桥蛋白轴 | TAM 分泌骨桥蛋白 SPP1，结合 CD44 促进肿瘤侵袭和基质重塑 |
| MIF–CD74 | 炎症调控 | 巨噬细胞迁移抑制因子，在 Visium 级别信号可检测 |

**揭示的生物学现象**：

该图将描述性的"Treg-TAM-CAF 共定位"结论升级为**机制性结论**——免疫抑制生态位中富集的信号通路揭示了细胞间通讯网络：CXCL12-CXCR4 和 CCL22-CCR4 是招募信号，TGF-β 和 IL-10 是维持信号，PD-1/TIGIT/LAG3 信号轴是效应 T 细胞耗竭的执行机制。

#### `niche_signature_volcano.png`

**内容**：火山图（Volcano Plot），横轴为 log2FC（niche_high vs niche_low，niche_high 为 q=0.85 截断），纵轴为 -log10(FDR)。
- 红色点：显著上调基因（log2FC > 0.5 且 FDR < 0.05）
- 蓝色点：显著下调基因
- 灰色点：不显著基因
- 选择性标注先验目标基因（FOXP3、TGFB1 等）

**如何理解**：右上角（高 log2FC、高 -log10FDR）的基因是在免疫抑制 niche 中最显著上调的候选分子。稀有免疫基因（如 FOXP3）由于 Visium 稀释效应，FC 值通常偏低，需配合 `niche_fraction_scatter.png` 一起解读。

#### `niche_fraction_scatter.png`

**内容**：检出率差值散点图，横轴为 log2FC，纵轴为 Δfrac（检出率差值 = frac_high - frac_low）。
- 红色点：双重显著（log2FC > 0.5 且 Δfrac > 0.10）
- 橙色点：仅检出率富集（Δfrac > 0.10 但 FC 不高）—— 典型稀有免疫基因（FOXP3）
- 蓝色点：仅 FC 高（可能是高表达管家基因）
- 灰色点：两者均不显著

**如何理解**：

本图以 Δfrac 为纵轴，直接展示"niche_high 相对 niche_low 有多少更多 spot 能检出该基因"，是火山图的重要补充视角。橙色区域（仅 Δfrac 显著）中的基因通常是低细胞丰度但具有功能重要性的稀有免疫细胞标志基因，在传统均值 log2FC 分析中往往被漏掉。

#### `sensitivity_niche_stability.png`

**内容**：敏感性分析结果图（v2 连续指标），展示两个参数维度（kNN 邻居数 / 半径倍增系数）的 Spearman ρ 和加权 Jaccard 两指标稳定性（2×2 子图），以 k=15 为参考配置。

**如何理解**：若各格子的数值接近 1（颜色深），则说明参数扰动对 niche 评分排序影响极小，结论稳健。阈值参考：Spearman ρ > 0.90 且 Weighted Jaccard > 0.80 为合格。

#### `param_scan_deg_stability_heatmap.png`

**内容**：k × niche_high_quantile 二维参数扫描热图，颜色编码不同参数组合下显著差异基因（FDR<0.05 且 log2FC>0.5）的数量。

**如何理解**：

- **横轴**：`niche_high_quantile`（0.70 / 0.75 / 0.80 / 0.85）
- **纵轴**：kNN 邻居数 k（8 / 10 / 15 / 20 / 25）
- **颜色**：颜色越深说明该参数组合下发现的显著 DEG 越多

选参标准：**颜色最深且处于"高原"区域**（与相邻格子结果相近）的参数组合为推荐。**当前主流程参数 k=15, q=0.85 位于高原区**，参数选择合理，结论稳健。

---

## 三、Step 3 — 论文图表生成（`pipeline/run_paper_figures.py`）

> 输出根目录：`results/paper_figures/`

Step 3 基于 Step 1（HCC4R AnnData）和 Step 2（niche scores CSV）的输出，生成适合期刊投稿的高质量图表（DPI=300），每张子图独立保存。通过 `python code/run_hcc4r.py --paper-figures` 或 `--figures-only` 触发。

> **修改绘图参数**：只需编辑 `code/pipeline/paper_plot_functions.py`，无需改动数据准备逻辑。

### Figure 1：整体课题设计与单细胞图谱

| 文件 | 内容 | 来源 |
|------|------|------|
| `fig1A_workflow_diagram.png` | 课题技术路线图（三列流程：scRNA-seq → 空间转录组 → 临床验证） | 自动生成（matplotlib 绘制，无数据文件依赖） |
| `fig1B_scrna_tsne_celltype.png` | scRNA-seq tSNE 细胞类型图 | 复制自 `results/HCC4R/scrna_tsne_celltype.png` |
| `fig1C_scrna_marker_heatmap.png` | 细胞类型 Marker 基因热图 | 复制自 `results/HCC4R/scrna_celltype_marker_heatmap.png` |
| `fig1D_t_cell_dotplot.png` | T/NK 子簇 Treg 标志基因气泡图 | 复制自 `results/HCC4R/t_cell_dotplot_horizontal.png` |

---

### Figure 2：Cell2location 空间反卷积结果

| 文件 | 内容 |
|------|------|
| `fig2A_spatial_celltype_all.png` | 所有主要细胞类型（Hepatocyte、Treg、T/NK、Myeloid、Fibroblast、HSC、B cell、Endothelial、Malignant）比例的 3×3 空间分布网格图，paper_ybp 配色 |
| `fig2B_spatial_hepatocyte.png` | Hepatocyte 比例空间分布（单图，论文主图质量） |
| `fig2C_spatial_treg.png` | Treg 比例空间分布（单图，论文主图质量） |
| `fig2D_spatial_myeloid.png` | Myeloid 比例空间分布（单图） |
| `fig2E_spatial_fibroblast.png` | Fibroblast 比例空间分布（单图） |

**如何理解**：这一组图将 Cell2location 反卷积的细胞类型估计结果可视化到二维组织切片坐标系上，直观展示了各细胞类型的空间分布格局。对照病理学切片图，可验证肿瘤核心（高 Hepatocyte）、免疫浸润区（高 T/NK、Myeloid）和基质区（高 Fibroblast）的物理位置与生物学预期是否一致。

---

### Figure 3：免疫抑制微生态位发现

| 文件 | 内容 |
|------|------|
| `fig3A_spatial_niche_clusters.png` | Leiden 邻域组成聚类结果空间分布（离散色板，每种 cluster 一色） |
| `fig3B_spatial_niche_semantic.png` | 语义 niche 标签图（紫色=`immunosuppressive_niche`，灰色=其他） |
| `fig3C_spatial_niche_score.png` | 综合免疫抑制 niche 评分空间热图（paper_ybp 渐变色，亮黄=高分） |
| `fig3D_spatial_niche_high.png` | niche_high（q=0.85 截断）二值分布图 |
| `fig3E_region_score_boxplots.png` | 四区域（tumor_core / tumor_edge / stroma_immune / other）免疫抑制评分箱线图 |
| `fig3F_hepatocyte_treg_scatter.png` | Hepatocyte vs Treg 比例散点图（颜色=niche评分，magma色标） |

**如何理解**：Figure 3 是本研究最核心的结论图组。从 3A（无监督发现）→ 3B（语义注释）→ 3C（连续评分）→ 3D（高分区定义）→ 3E/3F（区域比较与空间关系），完整展示了免疫抑制微生态位的识别过程和空间特征，逻辑自洽、层层递进。

---

### Figure 4：签名基因与细胞间通讯

| 文件 | 内容 |
|------|------|
| `fig4A_niche_volcano.png` | 签名基因火山图（log2FC vs -log10FDR，niche_high 基于 q=0.85 分组） |
| `fig4B_niche_fraction_scatter.png` | 检出率差值散点图（log2FC vs Δfrac，识别 FOXP3 等稀有免疫基因） |
| `fig4C_lr_communication_heatmap.png` | 配体-受体通讯分析热图（CXCL12-CXCR4、CCL22-CCR4 等关键信号轴） |
| `fig4D_celltype_correlation.png` | 细胞类型与 niche 评分相关性热图（vlag 配色，红=正相关，蓝=负相关） |

**如何理解**：Figure 4 将"描述性"的细胞空间分布结论升级为"机制性"的分子通讯证据。4A/4B 展示免疫抑制区域的差异表达基因，4C 揭示驱动 Treg 招募和维持的信号通路（招募-维持-耗竭三级网络），4D 量化细胞类型共定位的相关性。

---

### Figure 5：CHC20 验证队列（跨队列验证）

| 文件 | 内容 |
|------|------|
| `fig5A_chc20_spatial_niche_score.png` | CHC20 切片综合免疫抑制 niche 评分空间热图（paper_ybp 配色，独立验证） |
| `fig5B_chc20_spatial_niche_semantic.png` | CHC20 切片语义 niche 标签图（紫色=immunosuppressive_niche） |
| `fig5C_chc20_spatial_treg.png` | CHC20 切片 Treg 比例空间分布图 |
| `fig5D_rf_confusion_matrix.png` | Random Forest 标签迁移混淆矩阵（HCC4R niche labels → CHC20 spot 预测） |
| `fig5E_rf_feature_importance.png` | Random Forest 特征重要性柱状图（Top 15 细胞类型特征对 niche 识别的贡献） |

**如何理解**：

Figure 5 提供了免疫抑制微生态位的跨队列独立验证证据。CHC20 是完全独立的患者样本，在独立 Step 1+Step 2 分析后，若同样能识别出具有相似空间结构和细胞组成特征的免疫抑制生态位（5A/5B/5C），则证明该发现不依赖于特定样本。

Random Forest 验证（5D/5E）采用 HCC4R 数据训练模型，预测 CHC20 中各 spot 的 niche 标签：
- **混淆矩阵（5D）**：若 `immunosuppressive_niche` 的预测精度高（对角线值大），说明 CHC20 中存在与 HCC4R 相同特征的免疫抑制区域
- **特征重要性（5E）**：Treg、Myeloid 等细胞类型特征若位居前列，再次从机器学习角度证明这些细胞类型是免疫抑制生态位的核心定义要素

> CHC20 分析结果位于：`results/CHC20/spatial_niche/spatial_niche_scores.csv`

---

### Figure 6：参数稳健性验证

| 文件 | 内容 |
|------|------|
| `fig6A_sensitivity_stability.png` | kNN 邻居数和半径倍增系数敏感性分析（Spearman ρ + 加权 Jaccard，以 k=15 为参考） |
| `fig6B_param_scan_heatmap.png` | k × quantile 二维参数扫描热图（颜色=显著 DEG 数量，标注最优参数 k=15, q=0.85） |

**如何理解**：

Figure 6 是方法可靠性的自我验证图，回答"你的参数选择是任意的吗？"这一审稿人必问问题：
- **6A（敏感性分析）**：各参数扰动下 Spearman ρ > 0.90 且 Weighted Jaccard > 0.80，说明 niche 评分排序对参数选择不敏感
- **6B（参数扫描热图）**：k=15, q=0.85 位于颜色最深的"高原区"，说明该参数组合能最大化捕获显著 DEG，且不是孤立的局部最优点

---

## 四、Step 4（可选）— TCGA 生存分析（`pipeline/tcga_survival_analysis.R`）

基于 Step 2 生成的签名基因列表（`results/HCC4R/spatial_niche/immunosuppressive_niche_signature_genes.txt`），在 TCGA-LIHC 肝癌队列中进行 ssGSEA 评分和生存分析。本分析为可选步骤。

主要输出文件：

| 文件 | 内容 |
|------|------|
| `tcga_signature_score.csv` | 每位患者的 ssGSEA 免疫抑制评分 |
| `tcga_signature_survival.csv` | 评分 + 临床信息合并表 |
| `cox_results.txt` | Cox 比例风险回归结果（独立预后验证） |
| `km_plot.png` | Kaplan-Meier 生存曲线（高/低评分组对比） |
| `cox_forest_plot.png` | Cox 森林图（多变量独立性展示） |

**揭示的生物学现象**：

ssGSEA 分数是将空间发现的免疫抑制微生态位特征"投影"到大队列数据中的量化工具，相当于用一个精确的"免疫抑制程度仪"为每位肝癌患者打分。若高分组患者的生存时间显著短于低分组，则从流行病学角度证明了**肝癌空间免疫抑制微生态位的临床预后价值**——不仅是一种组织学现象，更是患者生死攸关的临床因素。

---

## 五、分析流程总结与各文件间的逻辑关系

```
data/scRNA_reference.h5ad
data/HCC4R/           (Visium Space Ranger 输出)  →  主分析 Discovery Cohort
data/CHC20_Visium/    (Visium Space Ranger 输出)  →  独立验证 Validation Cohort
         │
         ▼
code/run_hcc4r.py  →  pipeline/run_preprocessing.py (Step 1)
  ├── [results/HCC4R/]
  │     ├── adata_vis_post.h5ad                  ─── HCC4R 反卷积结果（核心 AnnData）
  │     ├── adata_sc_post.h5ad                   ─── scRNA 参考 RegressionModel 结果
  │     ├── spot_cell_proportion_HCC4R.csv        ─── HCC4R spot 细胞比例表
  │     ├── shared_genes_HCC4R.txt               ─── scRNA ∩ HCC4R 共有基因列表
  │     ├── scrna_tsne_celltype.png               ─── scRNA tSNE 细胞类型图 (Fig 1B)
  │     ├── scrna_celltype_marker_heatmap.png     ─── Marker 基因热图 (Fig 1C)
  │     ├── t_cell_dotplot_horizontal.png         ─── Treg 标志基因气泡图 (Fig 1D)
  │     └── regression_training_history_HCC4R.png ─── RegressionModel 训练曲线
  │
  │  （可选：同步运行 CHC20 Step 1）
  ├── [results/CHC20/]
  │     ├── adata_vis_post.h5ad                  ─── CHC20 独立反卷积结果
  │     ├── spot_cell_proportion_CHC20.csv
  │     ├── shared_genes_CHC20.txt
  │     └── regression_training_history_CHC20.png
  │
  ▼
pipeline/run_spatial_niche_analysis.py（Step 2，k=15，q=0.85）
  ├── [results/HCC4R/spatial_niche/]
  │     ├── spatial_niche_scores.csv                           ─── 核心评分表（含 niche_high 基于 q=0.85）
  │     ├── spatial_niche_parameters.csv                       ─── 分析参数记录（含 k=15, q=0.85）
  │     ├── sensitivity_analysis.csv                           ─── kNN/radius 稳健性（连续指标）
  │     ├── neighborhood_cluster_stats.csv                     ─── Leiden 聚类注释统计
  │     ├── immunosuppressive_niche_signature_genes_ranked.csv ─── 三层筛选签名基因排名表
  │     ├── immunosuppressive_niche_signature_genes.txt        ─── 签名基因列表（TCGA 接口）
  │     ├── prior_gene_set_auc.csv                             ─── Layer 2 先验基因 AUC 结果
  │     ├── gini_score_genes.csv                               ─── Layer 3 Gini 特异性基因（阈值0.3）
  │     ├── param_scan_deg_stability.csv                       ─── k×quantile 参数扫描数据（验证 k=15, q=0.85）
  │     └── plots/
  │           ├── spatial_hepatocyte.png                       ─── Hepatocyte 空间图
  │           ├── spatial_treg.png                             ─── Treg 空间图
  │           ├── spatial_myeloid.png                          ─── Myeloid 空间图
  │           ├── spatial_fibroblast.png                       ─── Fibroblast 空间图
  │           ├── spatial_immunosuppressive_niche_score.png    ─── 综合 niche 评分空间图
  │           ├── spatial_neighborhood_clusters.png            ─── Leiden 聚类空间图
  │           ├── spatial_niche_semantic_labels.png            ─── 语义 niche 标签图
  │           ├── spatial_region_labels.png                    ─── 辅助区域标注图
  │           ├── spatial_niche_high_score_spots.png           ─── niche_high(q=0.85) 二值图
  │           ├── spatial_niche_cluster_vs_score_comparison.png ─── 聚类法 vs 评分法对比
  │           ├── distance_to_hep_high_vs_niche_score.png      ─── 距离-评分折线图
  │           ├── region_score_boxplots.png                    ─── 各区域评分箱线图
  │           ├── hepatocyte_vs_treg_niche_score.png           ─── Hepatocyte-Treg 散点图
  │           ├── celltype_niche_correlation.png               ─── 细胞类型相关性热图
  │           ├── lr_communication_heatmap.png                 ─── L-R 通讯分析热图
  │           ├── niche_signature_volcano.png                  ─── 签名基因火山图
  │           ├── niche_fraction_scatter.png                   ─── 检出率差值散点图
  │           ├── sensitivity_niche_stability.png              ─── 稳健性分析图（Spearman ρ + Jaccard）
  │           └── param_scan_deg_stability_heatmap.png         ─── k×q 参数扫描热图（最优参数验证）
  │
  │  （同步：CHC20 独立 Step 2）
  └── [results/CHC20/spatial_niche/]   ─── 结构与 HCC4R 相同，用于 Figure 5 验证
  │
  ▼
pipeline/run_paper_figures.py（Step 3，可选）
  └── [results/paper_figures/]
        ├── fig1A_workflow_diagram.png
        ├── fig1B_scrna_tsne_celltype.png
        ├── fig1C_scrna_marker_heatmap.png
        ├── fig1D_t_cell_dotplot.png
        ├── fig2A_spatial_celltype_all.png
        ├── fig2B_spatial_hepatocyte.png  ─ fig2E_spatial_fibroblast.png
        ├── fig3A_spatial_niche_clusters.png ─ fig3F_hepatocyte_treg_scatter.png
        ├── fig4A_niche_volcano.png ─ fig4D_celltype_correlation.png
        ├── fig5A_chc20_spatial_niche_score.png ─ fig5E_rf_feature_importance.png
        └── fig6A_sensitivity_stability.png / fig6B_param_scan_heatmap.png
         │
         ▼
pipeline/tcga_survival_analysis.R（可选，Step 4）
  ├── tcga_signature_score.csv       ─── 每位患者的 ssGSEA 评分
  ├── tcga_signature_survival.csv    ─── 评分 + 临床信息合并表
  ├── cox_results.txt                ─── Cox 回归结果（独立预后验证）
  ├── km_plot.png                    ─── KM 生存曲线（直观预后差异）
  └── cox_forest_plot.png            ─── 森林图（多变量独立性展示）
```

**核心生物学结论链**：

1. **Cell2location 空间反卷积**（Step 1）：scRNA-seq 参考签名驱动的贝叶斯反卷积，将每个 Visium spot 的混合 RNA 信号分解为各细胞类型的绝对丰度估计，从"基因表达"转化为"细胞类型空间分布"
2. **空间 kNN 邻域构建（k=15，参数扫描最优值）**（Step 2 阶段1）：计算每个 spot 周围 15 个物理邻居的细胞组成均值，消除单点噪声，捕获真实的"局部微环境"特征
3. **Leiden 邻域聚类 + 多维语义注释**（Step 2 阶段2）：无监督发现细胞类型共定位的空间生态结构 → 基于 Treg/Myeloid/Fibroblast/基因评分四维综合排名，识别 `immunosuppressive_niche`
4. **连续评分 + 高分位数截断（q=0.85，参数扫描最优值）**（Step 2 附加）：为每个 spot 独立计算连续评分，以 85 百分位为阈值划分 `niche_high`/`niche_low`，用于签名基因的统计检验分组
5. **三层基因筛选策略**（Step 2）：Wilcoxon+FDR（Layer 1）+ 先验 AUC（Layer 2）+ Gini 指数（Layer 3，阈值0.3），弥补 Visium 稀释效应对稀有免疫基因（FOXP3 等）的不敏感问题，确保签名基因集的生物学完整性
6. **L-R 通讯分析**（Step 2）：揭示生态位内部的细胞间分子通讯网络（CXCL12-CXCR4 招募、TGFB1-TGFBR1 维持、PD-1/TIGIT/LAG3 耗竭），将共定位描述升级为机制性证据
7. **CHC20 独立跨队列验证 + Random Forest 标签迁移**（Step 3 Figure 5）：CHC20 独立分析证明 niche 的普遍性；RF 模型从 HCC4R 迁移至 CHC20 的高准确率证明两队列 niche 特征的一致性
8. **参数稳健性双重验证**（Step 3 Figure 6）：敏感性分析（Spearman ρ > 0.90）+ 参数扫描热图（k=15, q=0.85 位于高原区），从方法论层面保证结论不依赖参数的任意选择
9. **TCGA 大队列验证（可选）**（Step 4）：将空间转录组发现（两张切片）投影到 TCGA-LIHC 大队列（n≈370）→ ssGSEA 评分 + KM 曲线 + Cox 回归，证明免疫抑制微生态位具有普遍临床意义和独立预后价值
