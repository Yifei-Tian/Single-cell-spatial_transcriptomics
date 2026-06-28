# 需求文档

我计划用HCC4R数据作为主分析的空间转录组数据，CHC20作为验证分析的空间转录组数据，并做出一些能够插入论文中的图表。

要求：这里面的所有子图都要单独生成，不要多个子图拼接。

目前我的代码已经有这些图表的一部分，将尚未生成的图表生成。

## Figure 1: 课题总体设计与单细胞参考数据集的细胞图谱 (Baseline)

> **目的：** 阐明研究设计；利用 scRNA-seq 数据建立高质量的细胞类型参考，重点勾勒出免疫抑制细胞（如 Treg）的转录特征。

- ### **1A. 课题技术路线图 (Workflow Diagram)**

  - **内容：** 流程示意图。左边是 scRNA-seq（细胞分群与 Marker 鉴定），中间是 2 例空间转录组切片（HCC4R 作为发现集，CHC20 作为独立验证集，进行解卷积与空间生态位分析），右边是外部临床大队列（TCGA/ICGC 肝癌数据，进行生态位特征的生存预后验证）。

- ### **1B. scRNA-seq 细胞分群 tSNE / UMAP 图 (`scrna_tsne_celltype.png`)**

  - **内容：** 展示单细胞数据集中所有细胞的低维嵌入投影，按主要细胞类型（T/NK, Myeloid, B cell, Malignant, Endothelial, HSC 等）着色，中心有白色背景的类型标签。

- ### **1C. 细胞类型特异性 Marker 基因热图 (`scrna_celltype_marker_heatmap.png`)**

  - **内容：** Yellow-Black-Purple 配色。X 轴为细胞类型，Y 轴为特异性 Marker（如 Malignant 的 `ALB`/`APOA2`，HSC 的 `ACTA2`，T 细胞的 `CD3D` 等），展示清晰的对角线高亮模式，证明单细胞注释的绝对可靠。

- ### **1D. T/NK 细胞亚群细分及 Treg 鉴定气泡图 (DotPlot)**

  - **内容：** 聚焦 T 细胞内部，展示 `FOXP3`、`CTLA4`、`TIGIT`、`IKZF2` 等免疫抑制标志物在 Treg 亚群中的特异性高表达（气泡大小代表表达比例，颜色深浅代表表达量），为后续空间解卷积提供精准的“Treg 签名”。

## 🗺️ Figure 2: 空间转录组特征解卷积与多细胞成分的联合映射

> **目的：** 将单细胞定义的细胞类型概率准确投射到空间切片上，观察各类细胞在肝癌组织中的空间分布初貌。

- ### **2A. 空间切片 H&E 染色病理图与解剖标注**

  - **内容：** 发现集切片（HCC4R）的原始 H&E 图像，旁边附带病理医生或基于形态学标注的肿瘤核心、癌旁、边界等区域。

- ### **2B. 核心细胞类型空间丰度图 (Spatial Scatterpie / FeaturePlot)**

  - **内容：** 在切片二维空间坐标上，用颜色深浅展示恶性肝癌细胞（Malignant）、HSC、Treg、Myeloid 细胞的预测占比（Proportion）。

- ### **2C. 细胞空间共定位相关性热图 (Spatial Co-localization Heatmap)**

  - **内容：** 计算切片上任意两种细胞丰度在所有 Spot 间的 Pearson 相关系数。
  - **揭示现象：** 揭示哪些细胞倾向于“成双成对”地出现在同一个微环境中（例如：Treg 与 HSC、Malignant 呈现强正相关，暗示免疫抑制轴的物理临近）。

## 🏔️ Figure 3: 主分析切片（HCC4R）的空间生态位（Spatial Niche）聚类与功能定量

> **目的：** 引入空间结构信息（Step 6 Leiden 聚类），将 Spot 划分为不同的“空间生态位”（Spatial Domain），并计算高阶微环境评分，锁定“免疫抑制前线”。

- ### **3A. 基于空间邻域图的空间生态位聚类图 (Spatial Domain / Niche Map)**

  - **内容：** 展示 `HCC4R` 经过 kNN 空间邻接与 Leiden 聚类后的结果。切片被涂成 4-5 种颜色，对应不同的解剖功能区（如 Domain 1: Tumor Core, Domain 2: Tumor Edge, Domain 3: Stroma/Immune Hub）。

- ### **3B. 空间生态位细胞组成条形图 (Composition BarPlot)**

  - **内容：** X 轴为各个 Spatial Domain，Y 轴为 100% 堆叠条形图，展示每个生态位内部的细胞成分比例。
  - **揭示现象：** 明确展示 `Tumor Edge`（肿瘤边缘区）中同时富集了 Malignant、HSC 和 Treg 细胞，构成了一个“多细胞混战”的独特空间结构。

- ### **3C. 多层次微环境评分空间投影与小提琴图组合图 (Violin + Spatial Plot)**

  - **内容：** 展示 `Treg_like_score`、`immune_stroma_score` 以及核心的 **`immunosuppressive_niche_score`（免疫抑制生态位评分）** 在切片上的连续分布，以及在不同 Domain 间对比的小提琴图。
  - **揭示现象：** 定量证实免疫抑制评分在 `Tumor Edge` 达到峰值。

## 🧪 Figure 4: 免疫抑制生态位特征基因提取与机制挖掘

> **目的：** 通过你代码中的“三层筛选策略”（Step 14），找出到底是谁在驱动这个免疫抑制生态位，挖掘其背后的生物学机制。

- ### **4A. 空间差异基因火山图 / 气泡图 (Volcano Plot / DotPlot)**

  - **内容：** 对比 `Tumor Edge` 与其他生态位，展示显著上调的空间差异表达基因（SDEGs），高亮突显基质重塑（如 `COL1A1`）与免疫检查点相关基因。

- ### **4B. 空间生态位特异性特征基因热图 (Niche Marker Heatmap)**

  - **内容：** 展示筛选出的特征基因（Signature Genes）在不同空间 Domain 中的表达情况，确立代表该免疫抑制微环境的“基因集签名”。

- ### **4C. 空间生态位通路富集分析气泡图 (GSVA / GSEA Bubble Plot)**

  - **内容：** 展示各个空间生态位激活的经典通路。重点突出免疫抑制生态位（Tumor Edge）中 **TGF-beta 信号通路、EMT（上皮-间质转化）、Angiogenesis（血管生成）以及 IL6-JAK-STAT3** 通路的显著富集。

## 🛡️ Figure 5: 独立空间转录组切片（CHC20）的标签迁移与盲测验证

> **目的：** 实施你提出的“主分析+独立验证”策略，利用模型将 HCC4R 的发现推广到 CHC20，证明规律的泛化性。

- ### **5A. 验证集（CHC20）区域标签迁移预测图 (Domain Label Transfer Spatial Plot)**

  - **内容：** 展示利用 `HCC4R` 训练的随机森林分类器在 `CHC20` 切片上的预测结果。检查模型自动识别出的 `Tumor Edge` 是否在解剖学上依然位于肿瘤边界。

- ### **5B. 验证集微环境评分与特征基因的一致性对比图 (Validation Violin / Cross-Cohort Comparison)**

  - **内容：** 左右并排。展示 `HCC4R` 和 `CHC20` 两个切片在各自对应的生态位（如 Tumor Edge）中，核心免疫抑制特征基因和评分的表达趋势是否高度同步（$p$ 值的 Wilcoxon 检验）。

- ### **5C. 空间半径邻域敏感性分析与参数扫描图 (`sensitivity_analysis_grid.png`)**

  - **内容：** 对应你代码中的 Step 11 和 Step 15（k × quantile 网格搜索）。展示随着邻域半径（radius_multiplier）调整，该免疫抑制生态位的生物学信号依然稳定存在，证明结论不是算法调参的偶然产物。

## ⏳ Figure 6: 外部大临床队列的预后验证与临床意义 (Translational Value)

> **目的：** 照应标题中的“预后验证”，将空间生态位衍生出的特征基因集（Signature）带入百例级临床样本（如 TCGA-LIHC 或 ICGC 肝癌队列），证明该空间微环境不仅存在，而且深刻影响患者的生存。

- ### **6A. 临床队列生存分析生存曲线 (Kaplan-Meier Survival Curve)**

  - **内容：** 根据你在空间免疫抑制生态位中筛选出的特征基因集，对 TCGA/ICGC 患者进行单样本基因集富集分析（ssGSEA）打分，分为 High-score 和 Low-score 两组。展示两组的总体生存率（OS）或无进展生存率（PFS）差异。
  - **预期结果：** 免疫抑制生态位评分高的患者，生存预后显著极差（$p < 0.05$）。

- ### **6B. 独立预后因素的 Cox 回归分析森林图 (Forest Plot for Cox Regression)**

  - **内容：** 展示单因素和多因素 Cox 回归分析结果。将生态位特征评分与临床指标（Age, Gender, TNM Stage, Grade）进行协同回归。
  - **预期结果：** 证明该免疫抑制空间微环境评分是肝癌患者的**独立不良预后因素**（HR > 1）。

- ### **6C. 临床病理特征关联热图 / 箱线图 (Clinical Correlation Plot)**

  - **内容：** 展示该空间微环境评分在不同临床分期（Stage I-IV）或病理分级（Grade 1-4）中的表达差异。通常随着恶性程度增高，该微环境评分显著上升。