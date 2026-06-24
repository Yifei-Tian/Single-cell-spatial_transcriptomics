# 分析流程输出文件说明

本文档对本项目所有 `results/` 目录下生成的输出文件进行逐一解释，说明每个文件的内容格式、如何阅读和解读，以及其背后揭示的生物学意义。

分析流程共分五个步骤，输出文件按步骤组织。

---

## 一、Step 1 — 数据预处理与 Cell2location 反卷积（`run_preprocessing.py`）

### 1.1 `adata_vis_post.h5ad`

**文件格式**：HDF5 格式的 AnnData 对象（单细胞/空间组学标准存储格式）

**内容说明**：

这是 CHC20 切片经过 Cell2location 贝叶斯反卷积之后的核心数据文件，是后续所有空间分析的起点。文件内部包含：

- `adata.X`：每个 spot 的原始 count 表达矩阵（行 = spot，列 = 基因）
- `adata.obs`：每个 spot 的元数据，包含组织位置信息
- `adata.obsm["spatial"]`：每个 spot 的二维空间坐标（来自 Visium 芯片物理布局）
- `adata.obsm["means_cell_abundance_w_sf"]`：Cell2location 后验估计的细胞类型丰度矩阵（行 = spot，列 = 细胞类型），这是反卷积的主要输出

**如何理解**：

10x Visium 空间转录组每个 spot 直径约 55 μm，理论上可覆盖 1–10 个细胞，但实际上每个 spot 捕获的 RNA 来自多种细胞类型的混合信号。Cell2location 是一种贝叶斯层次模型，利用 scRNA-seq 参考数据中各细胞类型的基因表达特征（参考签名），通过负二项分布似然函数将每个 spot 的混合信号"反卷积"为各细胞类型的估计数量（后验均值 `means_cell_abundance_w_sf`）。

**揭示的生物学现象**：

反卷积结果直接回答了"每个空间位置存在哪些细胞类型"这一核心问题。通过观察 `means_cell_abundance_w_sf` 矩阵，可以看到：
- 肿瘤实质（Hepatocyte 丰度高）在组织切片中的分布
- 免疫细胞（Treg、T/NK、Myeloid）的浸润区域及其空间梯度
- 基质细胞（Fibroblast）形成的物理屏障位置

这是从"基因表达"到"细胞类型空间分布"的关键转化步骤。

---

### 1.2 `adata_vis_post_CHC23.h5ad`（等效：`adata_vis_chc23_post.h5ad`）

**文件格式**：同上，HDF5 格式的 AnnData 对象

**内容说明**：

CHC23 切片独立运行 Cell2location 的结果，结构与 `adata_vis_post.h5ad` 完全相同，但对应不同患者的肿瘤组织样本。

**如何理解与解读**：

CHC23 是独立验证切片，其 Cell2location 建模使用了从同一 scRNA-seq 数据集中独立训练的参考签名（而非直接借用 CHC20 的签名），保证了验证的统计独立性。

**揭示的生物学现象**：

对比 CHC20 和 CHC23 的反卷积结果，可以评估肝癌免疫微环境的**跨患者一致性**。若两个切片均显示 Treg 与 Myeloid 细胞在同一空间区域共富集，则支持免疫抑制微环境的普遍性，而非个例现象。

---

### 1.3 `spot_cell_proportion.csv`

**文件格式**：CSV 表格

**列说明**（每列代表一种细胞类型，每行代表一个 spot）：

| 列名 | 含义 |
|------|------|
| `spot_id` | Visium spot 的唯一标识符（如 `ACGCCTGACACGCGCT-1`） |
| `Hepatocyte` | 该 spot 中肝细胞比例（0–1） |
| `Treg` | 调节性 T 细胞比例 |
| `T/NK` | T 细胞和 NK 细胞的合并比例 |
| `Myeloid` | 髓系细胞（包括 TAM、MDSC 等）比例 |
| `Fibroblast` | 成纤维细胞（CAF）比例 |
| `B cell` | B 细胞比例 |
| ... | 其他细胞类型 |

所有细胞类型比例之和为 1（归一化后）。

**如何理解**：

这是将 `adata_vis_post.h5ad` 中的绝对丰度值归一化为相对比例的结果表。每行是一个 spot 的细胞组成"快照"，直接反映该空间位置的局部微环境细胞构成。

**揭示的生物学现象**：

- `Treg` 列的数值分布可揭示肿瘤组织中免疫抑制压力的空间分布格局
- `Hepatocyte` 高值区域对应肿瘤实质（肿瘤核心区），低值区域对应间质/免疫浸润区
- `Treg` 与 `Myeloid` 同时高值的 spot 提示可能的免疫抑制微生态位区域

---

### 1.4 `spot_cell_proportion_CHC23.csv`

**文件格式**：同 `spot_cell_proportion.csv`

与 CHC20 的比例表对应，用于跨切片对比分析。

---

### 1.5 `t_cell_dotplot_horizontal.png`

**文件格式**：PNG 图像

**内容说明**：

横版气泡图（Dot Plot），展示 T/NK 细胞各子簇（subclusters，分辨率 res.3 下的 Leiden 簇）中四个关键标志基因的表达情况：

- **CD3D**：所有 T 细胞的通用标志
- **CD4**：辅助性 T 细胞和 Treg 均表达的标志
- **FOXP3**：Treg 的核心转录因子，特异性最强
- **IL2RA（CD25）**：Treg 活化标志，在稳定 Treg 中高表达

气泡大小编码该基因在该簇中的**表达细胞比例**（percent expressed），气泡颜色深浅编码**平均表达量**（mean expression）。

**如何理解**：

横版布局将 x 轴设为聚类簇，y 轴设为基因，方便在基因数目少、簇数目多时并排比较。找到 FOXP3 和 IL2RA 同时呈深色大气泡的那一行对应的簇，即为 Treg 亚群。

**揭示的生物学现象**：

FOXP3 是 Treg 的"主调控转录因子"（master transcription factor），其特异性表达是鉴定 Treg 的金标准。图中清晰显示 FOXP3 仅在某一特定簇（如簇9）中高表达，表明该簇为真正的调节性 T 细胞亚群，而非普通辅助 T 细胞或效应 T 细胞，为后续空间 Treg 分析奠定细胞类型定义基础。

---

### 1.6 `regression_training_history_CHC20.png` / `regression_training_history_CHC23.png`

**文件格式**：PNG 图像

**内容说明**：

Cell2location 第一阶段 RegressionModel（单细胞参考签名学习）的训练损失曲线，x 轴为训练 epoch 数，y 轴为 ELBO（证据下界，Evidence Lower BOund）。

**如何理解**：

ELBO 是变分推断（Variational Inference）中的优化目标，其绝对值越大（负 ELBO 越小）说明模型对数据的拟合越好。曲线应在若干 epoch 后趋于平稳（收敛），若曲线仍在下降则说明训练轮次不足。

**揭示的生物学现象**：

收敛的训练曲线保证了细胞类型参考签名的可靠性——只有当 RegressionModel 充分学习到每种细胞类型的特异性基因表达模式后，下一阶段的空间反卷积才能准确区分不同细胞类型。

---

### 1.7 `spatial_mapping_training_history_CHC20.png` / `spatial_mapping_training_history_CHC23.png`

**文件格式**：PNG 图像

**内容说明**：

Cell2location 第二阶段（空间建模）的训练损失曲线，同样展示 ELBO 随训练轮次的变化。

**如何理解**：

与参考签名训练曲线类似，此图验证空间建模是否充分收敛。若曲线未收敛，spot 级别的细胞丰度估计将不可靠。

---

### 1.8 `cross_slice_comparison/` 目录（多切片一致性验证图）

#### `cross_slice_mean_proportion_comparison.png`

**内容**：CHC20 与 CHC23 各细胞类型**全切片平均比例**的并排条形图（蓝色=CHC20，橙色=CHC23）。

**如何理解**：若两个切片的细胞类型组成比例趋势相近（如 Hepatocyte 均为最高占比，Treg 均偏低但存在），则说明两个样本的整体微环境构成具有代表性，不存在严重的批次偏差。

**揭示的生物学现象**：肝癌肿瘤微环境的"细胞生态"在不同患者之间是否具有共同特征，即肝细胞癌 TME 免疫细胞组成的普遍性。

#### `cross_slice_treg_distribution.png`

**内容**：CHC20 与 CHC23 Treg 比例分布的核密度估计（KDE）对比曲线。

**如何理解**：x 轴为 Treg 比例，y 轴为密度。分布峰值、宽度和尾部形状相近则说明两切片的 Treg 浸润模式一致。

**揭示的生物学现象**：Treg 在肝癌组织中的浸润通常是稀疏但局灶性富集的（分布呈右尾），这种"局灶性免疫抑制"模式的跨样本一致性是免疫逃逸机制普遍存在的空间证据。

#### `cross_slice_celltype_boxplot.png`

**内容**：Treg、Myeloid、Fibroblast、Hepatocyte 四种关键细胞类型在两切片中比例分布的箱线图对比。

**如何理解**：箱体中位线、四分位距和异常值分布的相似性反映两切片细胞比例的整体一致性。

---

### 1.9 `shared_genes_CHC20.txt` / `shared_genes_CHC23.txt`

**文件格式**：纯文本，每行一个基因名

**内容说明**：scRNA-seq 参考数据与对应 Visium 切片之间共有的基因列表，即最终用于 Cell2location 建模的基因集。

**如何理解**：共有基因数量通常在 2,000–5,000 个之间（取决于 scRNA-seq 测序深度和 Visium 捕获效率）。过少的共有基因会降低反卷积的区分能力。

---

## 二、Step 2 — 空间免疫抑制生态位分析（`run_spatial_niche_analysis.py`）

### 2.1 `spatial_niche/spatial_niche_scores.csv`

**文件格式**：CSV 表格（行 = spot，列 = 多种评分和标签）

**关键列说明**：

| 列名 | 类型 | 含义 |
|------|------|------|
| `spot_id` | 字符串 | Spot 唯一标识符 |
| `spatial_x` / `spatial_y` | 浮点数 | 空间坐标（像素单位） |
| `Hepatocyte` / `Treg` / `Myeloid` / `Fibroblast` / `T/NK` | 浮点数（0–1） | 各细胞类型归一化比例 |
| `immunosuppressive_gene_score` | 浮点数 | 免疫抑制标志基因（FOXP3、IL2RA、CTLA4、TIGIT、LAG3、TGFB1、IL10 等）的 CP10K+log1p 归一化平均表达量 |
| `neighborhood_cluster` | 字符串（如"0","1",...） | Leiden 邻域组成聚类分配的 cluster 编号 |
| `niche_semantic_label` | 字符串 | 语义注释标签：`immunosuppressive_niche` 或 `cluster_X` |
| `Treg_like_score` | 浮点数 | Z-score(Treg比例) + Z-score(免疫抑制基因评分) 之和 |
| `immune_stroma_score` | 浮点数 | Treg + T/NK + Myeloid + Fibroblast 四项 Z-score 之和 |
| `immunosuppressive_niche_score` | 浮点数 | 综合免疫抑制 niche 评分（5项 Z-score 加权求和，含邻域 Hepatocyte 均值） |
| `niche_high` | 布尔值 | 该 spot 是否属于高免疫抑制区域（评分超过 80 百分位） |
| `spatial_region` | 字符串 | 规则化空间区域标注（`tumor_core`/`tumor_edge`/`stroma_immune`/`other`） |
| `hep_high` | 布尔值 | 该 spot 是否属于高 Hepatocyte 比例区域（超过 75 百分位） |
| `distance_to_hep_high` | 浮点数 | 到最近高肝细胞区域的欧式距离 |
| `neighbor_Hepatocyte` 等 | 浮点数 | 空间半径邻域内各细胞类型比例的均值（平滑后的局部微环境信息） |

**如何理解整个评分体系**：

本分析采用两阶段策略识别免疫抑制生态位：

**第一阶段（数据驱动发现）**：以每个 spot 的 k-NN 邻域（k=15）内各细胞类型比例均值作为特征向量，对这些"邻域组成向量"运行 Leiden 无监督聚类，识别出若干个具有相似局部细胞组成的 cellular neighborhoods（`neighborhood_cluster` 列）。

**第二阶段（功能评分注释）**：对每个 Leiden cluster 计算四个维度的均值（Treg比例、Myeloid比例、Fibroblast比例、免疫抑制基因评分），选取综合排名最高（最大免疫抑制特征）的 cluster 语义注释为 `immunosuppressive_niche`。

额外计算的综合评分（`immunosuppressive_niche_score`）用于后续签名基因提取时的 `niche_high`/`niche_low` 分组，与 niche 发现过程本身解耦。

**揭示的生物学现象**：

- `niche_semantic_label = "immunosuppressive_niche"` 的 spot 集合定义了**肿瘤免疫抑制微生态位**的空间范围，即 Treg、TAM（Myeloid）和 CAF（Fibroblast）三类细胞同时局部富集的区域
- 这类区域在生物学上对应肿瘤逃逸的"免疫沙漠"或"免疫排除"区域，是效应 T 细胞难以发挥杀伤功能的关键空间结构
- `spatial_region` 的分层（肿瘤核心→肿瘤边缘→免疫基质区域）反映了肿瘤组织内部微环境的空间异质性

---

### 2.2 `spatial_niche/spatial_niche_parameters.csv`

**文件格式**：CSV 表格（两列：`parameter` 和 `value`）

**内容说明**：记录本次分析使用的所有关键参数及运行时计算得到的统计量，包括：

| 参数 | 含义 |
|------|------|
| `n_neighbors_knn` | 邻域组成聚类时的 k-NN 邻居数（默认15） |
| `leiden_resolution` | Leiden 聚类分辨率（默认0.5） |
| `n_neighborhood_clusters` | 实际识别出的 Leiden 聚类数量 |
| `neighbor_radius_multiplier` | 空间半径邻域倍增系数（默认1.25） |
| `neighbor_radius` | 实际使用的搜索半径（单位：坐标像素） |
| `hep_high_threshold` | 高肝细胞区域的 Hepatocyte 比例阈值 |
| `niche_high_threshold` | niche_high 评分阈值（80百分位值） |
| `n_niche_high` | 被标记为 niche_high 的 spot 总数 |
| `marker_genes_used` | 实际用于计算免疫抑制基因评分的基因列表 |

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
| `n_niche_high` | 该参数下 niche_high spot 的数量 |
| `niche_pct` | niche_high 占全部 spot 的百分比 |
| `score_mean` / `score_std` | niche 评分的均值和标准差 |

**如何理解**：

敏感性分析扫描了 kNN 邻居数（k=10/15/20）和空间半径倍增系数（×1.0/1.25/1.5）共 6 种参数组合，记录每种设置下高免疫抑制区域的占比变化。

关键判断标准：若 `niche_pct` 在不同参数下的变化幅度小（如均在 15%–25% 范围内波动），则说明 niche 识别结果**对参数选择不敏感**，结论具有稳健性。

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
| `mean_high` | 该基因在 `niche_high` spot 中的平均 log1p 表达量 |
| `mean_low` | 该基因在 `niche_low` spot 中的平均 log1p 表达量 |
| `log2_fc` | log₂(mean_high+1) - log₂(mean_low+1)，衡量在免疫抑制区域的富集倍数 |

按 `log2_fc` 降序排列，默认保留前 80 个基因。

**如何理解**：

`log2_fc > 0` 的基因在免疫抑制生态位中特异性高表达；`log2_fc` 越大，该基因在免疫抑制区域的富集越强。这些基因构成了免疫抑制微生态位的**转录特征签名（Signature）**。

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

同时在 `results/spatial_signature_genes.txt` 保存一份同内容副本供 R 脚本读取。

---

### 2.7 `spatial_niche/plots/` 可视化图表

#### `spatial_hepatocyte.png`

**内容**：Hepatocyte 细胞比例在切片空间坐标中的分布，viridis 色标（深色=低，亮色=高）。

**如何理解**：高 Hepatocyte 比例区域（亮色）对应肿瘤实质，即肝癌细胞密集区；低比例区域对应肿瘤间质（免疫浸润区、基质区）。此图是定义肿瘤核心区域的空间依据。

#### `spatial_treg.png`

**内容**：Treg 细胞比例的空间分布图（viridis 色标）。

**如何理解**：Treg 在肝癌中的分布通常并非均匀弥散，而是局灶性富集在特定区域（如肿瘤边界或免疫浸润区）。识别 Treg 富集的"热点"区域是本研究的核心目标之一。

**揭示的生物学现象**：局灶性 Treg 富集区域通常对应肿瘤内部的免疫抑制微生态位，这些区域中效应 T 细胞的杀伤功能受到最强烈的抑制，是肿瘤免疫逃逸的关键空间节点。

#### `spatial_myeloid.png` / `spatial_fibroblast.png`

**内容**：Myeloid（髓系细胞）和 Fibroblast（成纤维细胞）比例的空间分布。

**揭示的生物学现象**：
- Myeloid 细胞（包括肿瘤相关巨噬细胞 TAM）通常与 Treg 共分布，构成 Treg 招募与维持的细胞环境（CCL22-CCR4 轴）
- Fibroblast 高比例区域对应肿瘤基质，CAF 通过分泌 CXCL12 等趋化因子和重塑细胞外基质参与免疫抑制

#### `spatial_immunosuppressive_niche_score.png`

**内容**：综合免疫抑制 niche 评分的空间分布（magma 色标，深紫色→亮黄色代表评分从低到高）。

**如何理解**：magma 色标的选择使高评分区域（免疫抑制热点）在视觉上非常突出。整个切片的"评分地图"直观揭示了免疫抑制压力的空间分布格局。

**揭示的生物学现象**：高评分区域即免疫抑制微生态位的空间轮廓，该区域通常位于肿瘤实质与免疫浸润区的交界地带（即免疫细胞进入肿瘤的"前哨区域"），是制定靶向免疫治疗策略时需要重点关注的空间区域。

#### `spatial_neighborhood_clusters.png`

**内容**：Leiden 邻域组成聚类结果的空间分布图，每个 cluster 用不同颜色标识（tab20 色板）。

**如何理解**：这是分析方法的"中间过程可视化"，展示无监督聚类对肿瘤组织空间结构的划分。空间上相邻且颜色相同的 spot 组成一个具有相似局部细胞组成的"邻域"。

**揭示的生物学现象**：聚类结果反映了组织的**空间层级结构**：不同 cluster 对应肿瘤的不同功能区域（如肿瘤核心、免疫浸润区、基质区等）。这种无监督发现的空间分区与病理学上的组织学分区往往高度一致，说明基因表达驱动的细胞组成分布确实遵循组织生理学规律。

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

#### `hepatocyte_vs_treg_niche_score.png`

**内容**：Hepatocyte 比例（x 轴）vs Treg 比例（y 轴）散点图，颜色编码综合 niche 评分（magma 色标），虚线标注 75 百分位分界线划分四象限。

**如何理解**：重点关注右上象限（Hepatocyte 高 + Treg 高）中颜色最深的点，这些 spot 是肿瘤实质与 Treg 共定位最明显的区域。

**揭示的生物学现象**：Treg 在肿瘤实质内部的富集（而非局限于免疫浸润区边缘）是肝癌免疫抑制的关键特征，意味着 Treg 已成功突破物理屏障进入肿瘤核心，对肿瘤内的效应 T 细胞形成直接抑制。

#### `distance_to_hep_high_vs_niche_score.png`

**内容**：距离高肝细胞区域（肿瘤核心）的距离（x 轴，等频分箱）vs 平均免疫抑制 niche 评分（y 轴）的折线图，附带标准误差棒。

**如何理解**：若折线呈现"距肿瘤核心越近，免疫抑制评分越高"的趋势，则说明免疫抑制信号具有**距离依赖性**（肿瘤-近端区域受到更强的免疫抑制压力）。反之若无明显趋势则提示免疫抑制呈弥散分布。

**揭示的生物学现象**：肿瘤实质边界区域的强免疫抑制信号与肿瘤招募 Treg 和分泌趋化因子（如 CXCL12）的机制一致，是肿瘤主动塑造免疫抑制微环境的空间证据。

#### `region_score_boxplots.png`

**内容**：四个空间区域（`tumor_core`、`tumor_edge`、`stroma_immune`、`other`）中 `Treg_like_score`（左图）和 `immunosuppressive_niche_score`（右图）的箱线图。

**如何理解**：若 `stroma_immune` 或 `tumor_edge` 区域的箱体中位线最高，说明免疫抑制核心集中在免疫浸润区而非纯肿瘤实质，提示 Treg 的主要作用位点在"肿瘤-免疫交界面"。

#### `celltype_niche_correlation.png`

**内容**：各细胞类型比例与 niche 评分之间的 Pearson 相关系数热图（vlag 配色：红色=正相关，蓝色=负相关）。

**如何理解**：关注 Treg 与 Myeloid 之间的相关系数、以及这两者与 Hepatocyte 之间的相关系数。

**揭示的生物学现象**：
- **Treg 与 Myeloid 正相关**：支持 TAM 通过分泌 CCL22 招募 Treg 的机制假设（TAM→Treg 信号轴）
- **Treg 与 Hepatocyte 正相关或负相关**：揭示肿瘤细胞与 Treg 的空间关联模式（肿瘤是否主动招募 Treg 入侵实质）
- **Fibroblast 与 Myeloid 正相关**：提示 CAF-TAM 共定位在基质区形成双重免疫屏障

#### `lr_communication_heatmap.png`

**内容**：配体-受体（L-R）通讯分析热图，共两个子图：

- **左图**：7 个关键 L-R 信号对（行）在 `Niche-High` vs `Niche-Low` 两组（列）中的归一化信号强度热图（YlOrRd 色标）
- **右图**：每个 L-R 对的 log₂FC（Niche-High/Niche-Low），红色=在免疫抑制 niche 中上调，蓝色=下调

**7 个关键 L-R 信号对及其生物学意义**：

| L-R 对 | 信号通路 | 生物学功能 |
|--------|----------|-----------|
| CXCL12–CXCR4 | SDF-1 轴 | TAM 和基质细胞分泌 CXCL12，招募表达 CXCR4 的 Treg 进入肿瘤，是肝癌 Treg 招募的核心机制 |
| CCL22–CCR4 | 趋化因子轴 | M2 巨噬细胞分泌 CCL22，结合 Treg 表面 CCR4 受体，实现定向招募 |
| TGFB1–TGFBR1 | TGF-β 轴 | TGF-β 是 Treg 产生和维持的关键细胞因子，同时直接抑制效应 T 细胞 |
| PD-1–PD-L1 | 免疫检查点 | 肿瘤细胞/Myeloid 细胞上的 PD-L1 结合 T 细胞表面 PD-1，直接触发效应 T 细胞功能失活 |
| IL10–IL10RA | 抗炎细胞因子轴 | Treg 和 TAM 分泌 IL-10，通过 IL10RA 信号在局部造成免疫耐受环境 |
| TIGIT–NECTIN2 | 免疫检查点 | Treg 和耗竭 T 细胞表达 TIGIT，与 APC 上的 NECTIN2 结合触发抑制信号 |
| LAG3–MHC-II | 免疫检查点 | LAG3 结合 MHC-II 分子，与 TIGIT/PD-1 协同参与 T 细胞耗竭 |

**如何理解**：

左图中颜色越深说明该 L-R 对的通讯信号（配体基因表达 × 受体基因表达的乘积均值）越强；右图红色柱子说明该信号在免疫抑制 niche 中显著富集（log₂FC > 0）。

**揭示的生物学现象**：

该图将描述性的"Treg-TAM-CAF 共定位"结论升级为**机制性结论**——免疫抑制生态位中富集的信号通路揭示了细胞间通讯网络：CXCL12-CXCR4 和 CCL22-CCR4 是招募信号，TGF-β 和 IL-10 是维持信号，PD-1/TIGIT/LAG3 信号轴是效应 T 细胞耗竭的执行机制。这一"招募-维持-耗竭"三级信号网络完整地刻画了肝癌免疫逃逸的分子机制。

#### `sensitivity_niche_stability.png`

**内容**：敏感性分析结果图（对应 `sensitivity_analysis.csv`），包含两个子图：左图为不同 kNN 邻居数下 niche-high spot 占比，右图为不同半径倍增系数下的占比，红色虚线标注主分析参数。

**如何理解**：若各条形图高度接近红色虚线且变化幅度小，则说明结论稳健。

---

## 三、Step 3 — 差异表达分析（`run_de_analysis.py`）

### 3.1 `spatial_signature_genes.txt`（也在 `spatial_niche/` 下有副本）

**文件格式**：纯文本，每行一个 HGNC 基因名，默认 50 个基因

**内容说明**：

使用 Wilcoxon 秩和检验，比较 `niche_high`（免疫抑制 spot）与 `niche_low`（其他 spot）之间的基因表达差异，提取在免疫抑制区域中特异高表达的前 50 个基因。

表达量矩阵经过 CP10K 归一化和 log1p 变换以消除测序深度偏差。Wilcoxon 检验是非参数检验，不假设数据服从正态分布，适用于单细胞/空间转录组数据的高噪声特性。

**如何理解**：

该文件是整个分析流程中连接**空间转录组发现**与**临床队列验证**的关键接口。基因列表直接输入 R 脚本用于 TCGA-LIHC 队列的 ssGSEA 打分。

**揭示的生物学现象**：

相比 Step 2 中基于 log2FC 的特征基因（`immunosuppressive_niche_signature_genes_ranked.csv`），本文件通过 Wilcoxon 统计检验进一步筛选，保留了在统计上显著高表达的基因，降低了噪声基因的干扰。这些基因代表了肝癌空间免疫抑制微生态位最可靠的分子标志物，是具有潜在临床转化价值的候选靶点或生物标志物。

---

## 四、Step 4–5 — TCGA 生存分析（`tcga_survival_analysis.R`）

### 4.1 `tcga_signature_score.csv`

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| 患者 ID（TCGA barcode） | 唯一患者标识符（如 `TCGA-BC-A10T`） |
| `score` | ssGSEA（单样本基因集富集分析）评分（通常在 -1 到 1 之间，经过归一化） |

**如何理解**：

ssGSEA 将空间转录组发现的签名基因集（50–80 个基因）作为"感兴趣的基因集"，对 TCGA-LIHC 队列中每位患者的 bulk RNA-seq 表达谱进行单样本富集分析。高 ssGSEA 评分意味着该患者肿瘤中签名基因集整体表达水平高，即其肿瘤微环境具有**与空间转录组发现的免疫抑制生态位相似的分子特征**。

**揭示的生物学现象**：

ssGSEA 分数是将空间发现的免疫抑制微生态位特征"投影"到大队列数据中的量化工具，相当于用一个精确的"免疫抑制程度仪"为每位肝癌患者打分，分数越高代表其肿瘤承受的免疫抑制压力越大。

---

### 4.2 `tcga_signature_survival.csv`

**文件格式**：CSV 表格

**列说明**：

| 列名 | 含义 |
|------|------|
| 患者 ID | TCGA barcode |
| `score` | ssGSEA 免疫抑制评分 |
| `OS` | 总生存时间（overall survival，单位：天） |
| `OS.event` | 生存状态（1=死亡，0=截尾/存活） |
| `score_group` | 按中位数或最优截点分组（`high`/`low`） |
| `age` / `stage` 等 | 临床协变量（用于 Cox 多变量回归校正） |

**如何理解**：

该文件合并了 ssGSEA 评分与 TCGA 临床随访数据，是生存分析的输入数据表。

**揭示的生物学现象**：

通过分组变量 `score_group`，可以比较"免疫抑制高评分组"与"低评分组"患者的生存差异，直接验证"肝癌免疫抑制微环境越强，患者预后越差"这一核心假设。

---

### 4.3 `cox_results.txt`

**文件格式**：纯文本，R `summary(coxph(...))` 的标准输出格式

**内容说明**：

多变量 Cox 比例风险模型的回归结果摘要，包含：
- 各协变量（ssGSEA score、年龄、肿瘤分期）的回归系数（`coef`）
- 风险比（`exp(coef)`，即 Hazard Ratio，HR）
- 95% 置信区间
- Wald 统计量和 p 值
- 全局模型检验（对数秩检验、似然比检验、Wald 检验）

**如何理解**：

重点关注 ssGSEA score 对应行的 `exp(coef)`（HR）和 p 值：

- **HR > 1 且 p < 0.05**：免疫抑制评分每增加一个单位，死亡风险增加 (HR-1)×100%，该结论在排除年龄和分期等混杂因素后仍成立，说明免疫抑制微生态位特征是**独立预后因子**
- **HR > 1 但 p > 0.05**：信号方向正确但统计功效不足，可能需要更大样本
- **HR ≈ 1**：特征基因集与预后无显著关联，需要重新审查签名基因的选择

**揭示的生物学现象**：

多变量 Cox 模型的核心价值在于**独立性检验**。仅有单变量分析（KM曲线）的结论可能被年龄、分期等混杂因素解释，而 Cox 多变量分析控制这些协变量后仍显著，则证明免疫抑制 niche 特征是真正独立的预后驱动因素，具有潜在的**临床转化价值**（如作为预测免疫治疗应答率的生物标志物）。

---

### 4.4 `km_plot.png`

**文件格式**：PNG 图像

**内容说明**：

Kaplan-Meier 生存曲线，将 TCGA-LIHC 患者按 ssGSEA 评分分为高分组（`score_high`，红色）和低分组（`score_low`，蓝色），绘制两组的累积生存概率随时间的变化曲线。图中通常包含：
- 对数秩（log-rank）检验 p 值
- 风险集（at risk）表格
- 中位生存时间

**如何理解**：

若高分组曲线始终低于低分组曲线（即高分组生存率更低），且 p < 0.05，则视觉上直接证明免疫抑制程度高与更差的预后相关。曲线间隔距离越大，预后差异越显著。

**揭示的生物学现象**：

KM 曲线是临床研究中最直观的生存差异展示方式。若高免疫抑制评分组的中位生存时间显著短于低分组，则从流行病学角度证明了**肝癌空间免疫抑制微生态位的临床预后价值**——不仅是一种组织学现象，更是患者生死攸关的临床因素。

这一结论的转化意义在于：未来可能通过检测肿瘤活检的免疫抑制相关基因表达谱（如通过 RT-qPCR 或 NanoString 技术），预测患者对免疫检查点抑制剂（PD-1/PD-L1 抗体）的应答概率，指导个性化治疗决策。

---

### 4.5 `cox_forest_plot.png`

**文件格式**：PNG 图像

**内容说明**：

Cox 回归森林图（Forest Plot），以图形化方式展示 `cox_results.txt` 中各协变量的风险比：
- 每行代表一个协变量
- 矩形（方块）代表点估计 HR，大小通常与样本量/权重成正比
- 横线代表 95% 置信区间
- 竖虚线在 HR=1 处标注（无效假设）
- 若置信区间整体在 HR=1 右侧（不跨越 1），则该因素的效应有统计学意义

**如何理解**：

快速查看：若 ssGSEA score 对应的置信区间水平线完全位于 HR=1 的竖线右侧，则直接在图上确认其统计显著性。图中同时展示年龄、分期的效应，可直观比较各因素的独立效应大小。

**揭示的生物学现象**：

若森林图显示 ssGSEA score 的 HR 和统计显著性与年龄、分期相当或更强，则有力支持免疫抑制微生态位特征作为肝癌独立预后因子的临床价值。

---

## 五、分析流程总结与各文件间的逻辑关系

```
pre.py
  └── data/scRNA_reference.h5ad
  └── data/chc20_visium.h5ad / chc23_visium.h5ad
          │
          ▼
run_preprocessing.py（Step 1）
  ├── adata_vis_post.h5ad            ─── 反卷积结果（CHC20主分析）
  ├── adata_vis_post_CHC23.h5ad      ─── 反卷积结果（CHC23验证）
  ├── spot_cell_proportion.csv       ─── 细胞类型比例表（CHC20）
  ├── spot_cell_proportion_CHC23.csv ─── 细胞类型比例表（CHC23）
  ├── t_cell_dotplot_horizontal.png  ─── Treg 鉴定图
  ├── cross_slice_comparison/        ─── 多切片一致性验证图组
  └── [训练曲线图]
          │
          ▼
run_spatial_niche_analysis.py（Step 2）
  ├── spatial_niche/spatial_niche_scores.csv           ─── 核心评分表
  ├── spatial_niche/neighborhood_cluster_stats.csv     ─── Leiden 聚类注释统计
  ├── spatial_niche/spatial_niche_parameters.csv       ─── 分析参数记录
  ├── spatial_niche/sensitivity_analysis.csv           ─── 稳健性验证
  ├── spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv ─── 签名基因（带排名）
  ├── spatial_niche/immunosuppressive_niche_signature_genes.txt        ─── 签名基因列表
  ├── spatial_signature_genes.txt                      ─── 副本（供 R 脚本读取）
  └── spatial_niche/plots/                             ─── 全套空间可视化图
          │
          ▼
colocation.py + run_de_analysis.py（Step 3）
  ├── spot_with_coloc_label.csv      ─── 共定位标签（供 DE 分析输入）
  └── spatial_signature_genes.txt    ─── 统计筛选后的签名基因（更新/确认）
          │
          ▼
tcga_survival_analysis.R（Step 4–5）
  ├── tcga_signature_score.csv       ─── 每位患者的 ssGSEA 评分
  ├── tcga_signature_survival.csv    ─── 评分 + 临床信息合并表
  ├── cox_results.txt                ─── Cox 回归结果（独立预后验证）
  ├── km_plot.png                    ─── KM 生存曲线（直观预后差异）
  └── cox_forest_plot.png            ─── 森林图（多变量独立性展示）
```

**核心生物学结论链**：

1. **Cell2location 反卷积**：将混合信号的 Visium spot 分解为细胞类型构成 → 实现空间分辨率的细胞类型定位
2. **邻域组成聚类**：发现细胞类型共定位的空间生态结构 → 识别 Treg-TAM-CAF 共富集的免疫抑制生态位
3. **L-R 通讯分析**：揭示生态位内部的细胞间分子通讯 → 将共定位描述升级为 CXCL12-CXCR4、CCL22-CCR4 等具体信号通路的机制证据
4. **特征基因提取**：将空间发现转化为可量化的基因签名 → 为大队列验证提供分子接口
5. **TCGA 验证**：将空间转录组（n=2 切片）发现投影到大队列（n≈370） → 证明免疫抑制微生态位具有普遍临床意义和独立预后价值
