# 新增结果文件详细说明

> 本文档对应本轮代码改进（2026-06 版本）所新增/修改的全部输出文件。
> 原有文件说明请参见 `docs/result_explain.md`。

---

## 目录

1. [run_preprocessing.py 新增改动](#1-run_preprocessingpy-新增改动)
2. [run_spatial_niche_analysis.py 新增结果](#2-run_spatial_niche_analysispy-新增结果)
   - 2.1 [敏感性分析（改进版）](#21-敏感性分析改进版)
   - 2.2 [特征基因提取（三层升级版）](#22-特征基因提取三层升级版)
   - 2.3 [L-R 通讯分析（改进版）](#23-l-r-通讯分析改进版)
   - 2.4 [niche 空间二值分布图（新增）](#24-niche-空间二值分布图新增)
   - 2.5 [聚类 vs 评分并排对比图（新增）](#25-聚类-vs-评分并排对比图新增)
   - 2.6 [参数扫描热图（新增 Step 15）](#26-参数扫描热图新增-step-15)
3. [run_chc23_validation.py 所有输出](#3-run_chc23_validationpy-所有输出)
   - 3.1 [CHC23 平行 Niche 分析结果](#31-chc23-平行-niche-分析结果)
   - 3.2 [跨切片一致性量化结果](#32-跨切片一致性量化结果)
   - 3.3 [CHC23 DE 验证结果](#33-chc23-de-验证结果)

---

## 1. run_preprocessing.py 新增改动

### 改动说明

在 `main()` 的 Step 7（CHC23 单细胞模型训练）中，CHC23 现在使用更保守的早停参数：

| 参数 | 原值 | 新值 | 变化原因 |
|------|------|------|----------|
| `max_epochs` | 250 | **400** | 给模型足够的训练时间 |
| `early_stopping_patience` | 30 | **50** | 扩大早停观察窗口，避免在短暂 ELBO 平台期误判收敛 |
| `early_stopping_min_delta` | 1e-4 | **5e-5** | 对微小改善更敏感，防止过早停止 |

### 影响的输出文件

| 文件 | 路径 | 说明 |
|------|------|------|
| CHC23 训练曲线 | `results/regression_training_history_CHC23.png` | 若训练更充分，曲线将呈现更长的下降段再趋于平稳 |
| CHC23 空间数据 | `results/adata_vis_post_CHC23.h5ad` | Cell2location 后验丰度分布应更具区分度（箱线图不再退化为一条线） |

### 如何判断是否修复成功

重新运行 `run_preprocessing.py` 后，查看 `regression_training_history_CHC23.png`：
- **成功**：ELBO 曲线在更晚的 epoch（>100 轮）才趋于平稳，最终值低于旧版本
- **成功**：`adata_vis_post_CHC23.h5ad` 中各细胞类型的后验丰度箱线图呈现明显差异
- **失败**：若 ELBO 仍在 ~90 轮停止，考虑进一步增大 `max_epochs` 至 600 或手动禁用早停

---

## 2. run_spatial_niche_analysis.py 新增结果

### 2.1 敏感性分析（改进版）

**文件路径：** `results/spatial_niche/sensitivity_analysis.csv`  
**对应图片：** `results/spatial_niche/plots/sensitivity_niche_stability.png`

#### 问题背景

原版本敏感性分析记录各参数设置下的 `niche_pct`（niche_high spot 占比），但由于使用固定分位数阈值（如 80%）切割 niche_high，导致 `niche_pct` 在所有参数下**永远等于 20%**，图表完全无法区分参数优劣。

#### 改进方案

新版本改用 **Jaccard 相似度** 作为稳定性指标：

- 以主分析参数（k=15 kNN）产生的 niche_high spot 集合作为**参考集合**
- 对每种其他参数设置，计算其选出的 niche_high 集合与参考集合的 Jaccard 系数：
  ```
  Jaccard = |A ∩ B| / |A ∪ B|
  ```

#### CSV 文件列说明

| 列名 | 含义 | 典型值 |
|------|------|--------|
| `param_mode` | 参数类型（`knn` 或 `radius`） | `knn` |
| `param_label` | 参数描述 | `kNN k=10` |
| `is_reference` | 是否为参考参数（k=15） | `True`/`False` |
| `n_niche_high` | niche_high spot 数量 | ~500 |
| `niche_pct` | niche_high spot 占比（%) | ~20.0 |
| `jaccard_vs_ref` | 与参考集合（k=15）的 Jaccard 相似度 | 0.72–1.00 |
| `score_std` | niche 评分的标准差 | 任意正数 |

#### 结果解读

| Jaccard 值 | 含义 |
|-----------|------|
| > 0.85 | 参数选出的空间区域与主分析高度一致，**结果稳健** |
| 0.70–0.85 | 有一定差异，可能存在参数依赖性，建议关注 |
| < 0.70 | 参数变化对结果影响较大，需重新审视参数选择 |

#### 可视化图说明

`sensitivity_niche_stability.png`：
- 左图：k-NN 邻居数（10/15/20）vs Jaccard 相似度（红色柱 = 参考参数 k=15）
- 右图：半径倍增系数（×1.0/1.25/1.5）vs Jaccard 相似度
- 橙色虚线：Jaccard=0.85 的稳健性标准线
- 所有柱均接近 1.0 说明结果对邻域参数**不敏感**，是理想结果

---

### 2.2 特征基因提取（三层升级版）

原版本仅按 log2FC 排序，新版本引入三层并行筛选策略。

---

#### Layer 1：Wilcoxon + FDR 主签名基因

**文件路径：** `results/spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv`（格式升级）

##### 升级内容

原文件仅有 `gene / mean_high / mean_low / log2_fc` 四列，新版本新增：

| 新增列 | 含义 |
|--------|------|
| `pvalue` | Mann-Whitney U 检验 p 值（单侧，niche_high > niche_low） |
| `fdr` | Benjamini-Hochberg 方法校正后的 FDR |

##### 筛选逻辑变化

| 筛选维度 | 原版本 | 新版本 |
|----------|--------|--------|
| 统计检验 | 无 | Wilcoxon 秩和检验 |
| 主筛选条件 | 仅 log2FC > 0 | FDR < 0.05 且 log2FC > 0.5 |
| 补充条件 | 无 | 若严格条件不足 Top-N，放宽至 FDR < 0.2 且 log2FC > 0.3 |

##### 结果解读

- **FDR < 0.05 的基因**：在 niche_high vs niche_low 中表达差异具有统计学显著性，可安全用于报告
- **FDR 0.05–0.2 的基因**：有一定差异趋势，作为候选基因，需独立验证
- **log2FC 含义**：log2FC = 1 意味着 niche_high 的均值表达量是 niche_low 的 2 倍

---

#### Layer 1 衍生：火山图

**文件路径：** `results/spatial_niche/plots/niche_signature_volcano.png`

##### 图形说明

- **横轴**：log₂FC（niche_high / niche_low），正值表示 niche_high 中更高表达
- **纵轴**：-log₁₀(FDR)，越大表示越显著
- **颜色**：
  - 🔴 红色：显著上调（log2FC > 0.5 且 FDR < 0.05）
  - 🔵 蓝色：显著下调（log2FC < -0.5 且 FDR < 0.05）
  - ⚫ 灰色：不显著
- **标注**：FOXP3、TGFB1、FAP、CCL22、SPP1 等目标免疫抑制基因会自动标注名称
- **虚线**：垂直线 = |log2FC| = 0.5 的阈值；水平线 = FDR = 0.05 的阈值

##### 如何使用

若目标基因（如 FOXP3）出现在右上象限（红色且有标注），说明在空间上有显著的 niche 特异性表达，适合作为免疫抑制 niche 的生物标志物。

---

#### Layer 2：先验功能基因集 AUC 检验

**文件路径：** `results/spatial_niche/prior_gene_set_auc.csv`

##### 背景

纯 log2FC 筛选的缺陷在于：低表达稀有基因（如 FOXP3）因 Visium spot 级别的细胞稀释效应，表达量极低、log2FC 接近 0，永远无法进入 Top-50，但其在少数 niche_high spot 中的**检出率**可能显著高于 niche_low。AUC 能捕捉这种差异。

##### 检验内容

对三组先验功能基因集分别报告每个基因的统计指标：

| 基因集 | 代表基因 | 生物学意义 |
|--------|----------|-----------|
| Treg_markers | FOXP3, IL2RA, CTLA4, TIGIT, IKZF2 | Treg 细胞核心标志 |
| TAM_features | CD163, MRC1, TGFB1, IL10, CXCL12 | 肿瘤相关巨噬细胞特征 |
| CAF_activation | FAP, ACTA2, POSTN, COL1A1, CCL22 | 肿瘤相关成纤维细胞活化 |

##### CSV 列说明

| 列名 | 含义 |
|------|------|
| `gene_set` | 所属功能基因集（Treg/TAM/CAF） |
| `gene` | 目标基因名 |
| `actual_gene` | 实际使用的基因名（若主基因缺失则为替代基因）|
| `mean_high` | niche_high spot 中的均值表达量 |
| `mean_low` | niche_low spot 中的均值表达量 |
| `log2_fc` | log2 倍数变化 |
| `pvalue` | Mann-Whitney U 检验 p 值 |
| `fdr` | BH 校正后 FDR |
| `auc` | AUC（0.5 = 无区分力；> 0.6 = 有意义；> 0.7 = 良好） |

##### 结果解读

- **AUC > 0.6 且 FDR < 0.05**：该基因在 niche_high spot 中显著富集，可作为核心 signature 基因
- **AUC 0.5–0.6 且 FDR < 0.05**：有统计显著差异但区分力有限，作为辅助支持证据
- **AUC > 0.6 但 FDR > 0.05**：样本量可能不足或表达高度变异，谨慎使用

##### 如果 FOXP3 未达标

若 FOXP3 的 AUC < 0.6 或 FDR > 0.05，生物学解释为：
> FOXP3 在 Visium spot 级别高度稀释（每个 spot 包含数百个细胞），即使 Treg 富集，其 FOXP3 信号也被非 Treg 细胞的背景稀释。此时需借助 TIGIT、IKZF2 等在 Treg 以外也有一定表达的基因作为代理标志，或使用 cell2location 丰度值代替 FOXP3 转录信号。

---

#### Layer 3：Gini Index 特异性评分

**文件路径：** `results/spatial_niche/gini_score_genes.csv`

##### 原理

Gini Index（基尼系数）衡量数值分布的不均匀性：
- Gini = 0：所有 spot 中表达量完全均等
- Gini = 1：仅一个 spot 中有表达量，其余为 0

**高 Gini + niche_high 中有局灶性高表达 = 免疫抑制 niche 中高度特异的稀有基因**（如 FOXP3 集中在 Treg 密集区域的少数 spot）

##### CSV 列说明

| 列名 | 含义 |
|------|------|
| `gene` | 基因名 |
| `gini` | Gini 系数（在 niche_high spot 中计算） |
| `log2_fc` | log2 倍数变化（仅保留 log2FC > 0 的基因） |
| `mean_high` | niche_high spot 中的均值表达量 |

- 文件按 `gini` 降序排列，最顶部为最局灶性表达的基因
- 仅包含 `gini > 0.5` 且 `log2FC > 0` 的基因（避免非特异性高表达基因）

##### 三层方案联合解读建议

| 分类 | 条件 |
|------|------|
| **核心 signature 基因** | 三层均命中（FDR<0.05 + AUC>0.6 + Gini>0.5） |
| **强候选基因** | 任意两层命中 |
| **弱候选基因** | 仅一层命中，需独立验证 |

---

### 2.3 L-R 通讯分析（改进版）

**文件路径：** `results/spatial_niche/plots/lr_communication_heatmap.png`

#### 改进内容

##### 问题

原版本直接计算 `product = ligand × receptor`，因任一基因在 Visium spot 中接近 0（细胞稀释效应），乘积同样接近 0，导致 CCL22-CCR4、IL10-IL10RA 等 Treg 特异通路整行为 0，热图失去对比性。

##### 改进方案：秩归一化乘积（Rank-Normalized Product）

新版本对每个基因的表达向量做**秩归一化**（rankdata / n_spots），将值域映射到 [0, 1]，再计算配体秩 × 受体秩：
```
lig_ranked = rankdata(ligand_expression) / n_spots
rec_ranked = rankdata(receptor_expression) / n_spots
product = lig_ranked × rec_ranked  ∈ [0, 1]
```

优点：
- 解决稀疏零值导致乘积为 0 的问题
- 保留相对表达高低的顺序信息
- 兼容 log2FC 计算

##### 新增 L-R 通路

除原有 7 个通路外，新增以下 Visium 级别可检测的高表达通路：

| 通路 | 生物学意义 | 检测优势 |
|------|-----------|----------|
| SPP1–CD44 | TAM 分泌骨桥蛋白，招募并维持免疫抑制微环境 | SPP1 在 TAM 中高表达，信号强 |
| MIF–CD74 | 巨噬细胞迁移抑制因子，促进 TAM 极化 | MIF 广泛表达，检出率高 |
| VEGFA–KDR | 血管生成，与 CAF 活化相关 | VEGFA 在肿瘤基质中高表达 |
| Galectin9–TIM-3 | 抑制性检查点，LGALS9 诱导 T 细胞耗竭 | Visium 级别可检测 |
| CCL2–CCR2 | 单核/巨噬细胞招募 | CCL2 在 CAF 中高表达 |

##### GENE_FALLBACKS 替代基因机制

当主基因在数据集中缺失时，自动尝试同家族替代基因，并在热图标签中用 `*` 标注：
- CCR4 缺失 → 尝试 CCR2, CCR5
- IL10RA 缺失 → 尝试 IL10RB
- HAVCR2 缺失 → 尝试 TIM3, TIMD4

##### 图形解读

- 左图热图：行为 L-R 对，列为 Niche-High/Niche-Low，颜色编码秩归一化乘积均值（越亮表示通讯信号越强）
- 右图条形图：log2FC（Niche-High / Niche-Low），红色 = niche_high 中信号更强（上调），蓝色 = 下调
- `*` 标注的行使用了替代基因

---

### 2.4 niche 空间二值分布图（新增）

**文件路径：** `results/spatial_niche/plots/spatial_niche_high_score_spots.png`

#### 内容说明

基于综合免疫抑制 niche 评分（`immunosuppressive_niche_score`）的分位数阈值（默认 80th percentile）将所有 spot 切割为两类，并在空间坐标系上展示：

| 颜色 | 含义 |
|------|------|
| 🔴 红色 | 评分高于 80th percentile（`niche_high = True`），即综合评分最高的 20% spot |
| ⚫ 灰色 | 其余 spot（`niche_high = False`） |

#### 与聚类图的区别

此图与 `spatial_niche_semantic_labels.png`（Leiden 聚类方法）的核心区别：
- **聚类图**：基于每个 spot 的**邻域组成特征**（周围 k=15 个 spot 中的细胞类型比例）进行聚类
- **评分图**（本图）：基于每个 spot **自身的综合评分**（Treg 比例 + Myeloid 比例 + Fibroblast 比例 + 基因评分）进行阈值切割

> 两图之间的分歧有生物学意义——详见下一节。

---

### 2.5 聚类 vs 评分并排对比图（新增）

**文件路径：** `results/spatial_niche/plots/spatial_niche_cluster_vs_score_comparison.png`

#### 内容说明

同一画面内并排展示两种 niche 识别方法的空间分布，图题中自动标注两者的 Jaccard 相似度。

| 图 | 方法 | 数据来源 |
|---|------|---------|
| 左图 | Leiden 聚类注释 | `niche_semantic_label` 列中含 "immunosuppressive" 的 spot |
| 右图 | 评分阈值切割 | `niche_high = True` 的 spot |

#### 结果解读

| 情形 | 描述 | 生物学含义 |
|------|------|-----------|
| **高度重叠（Jaccard > 0.7）** | 两张图中红色区域基本一致 | 两种方法均指向同一空间区域，结论高度可靠，可作为论文核心证据 |
| **聚类有但评分无的 spot（边缘区域）** | 左图红色但右图灰色 | 这些 spot 的**邻域组成**像免疫抑制 niche（周围有 Treg/TAM/CAF），但**自身评分**（基因表达 + 细胞比例）不够高 → 免疫抑制生态位的"社会性"成分（受邻居影响而不是自身驱动） |
| **评分有但聚类无的 spot（孤立岛）** | 右图红色但左图灰色 | 这些 spot 的自身综合评分很高，但**周围邻域细胞组成**不典型 → 孤立的高免疫抑制 spot，生态位的"自主性"成分（不依赖邻域环境的局灶性免疫抑制） |
| **两者均无的 spot** | 两图均为灰色 | 普通肿瘤实质或基质区域，免疫抑制信号弱 |

---

### 2.6 参数扫描热图（新增 Step 15）

**CSV 文件：** `results/spatial_niche/param_scan_deg_stability.csv`  
**热图文件：** `results/spatial_niche/plots/param_scan_deg_stability_heatmap.png`

#### 背景

Leiden 分辨率（resolution）和 niche_high 分位数阈值（quantile）对分析结果有重要影响，但如何选择最优参数组合缺乏客观依据。参数扫描通过遍历参数网格，以 DEG 显著性作为目标函数，寻找最优参数区域。

#### 扫描范围

| 参数 | 扫描范围 | 步长 |
|------|---------|------|
| `leiden_resolution` | 0.3, 0.4, 0.5, 0.6, 0.7 | 0.1 |
| `niche_high_quantile` | 0.70, 0.75, 0.80, 0.85 | 0.05 |
| **总组合数** | **20 种** | - |

#### CSV 列说明

| 列名 | 含义 |
|------|------|
| `resolution` | Leiden 聚类分辨率 |
| `quantile` | niche_high 分位数阈值 |
| `n_clusters` | 该参数下 Leiden 识别的邻域聚类数量 |
| `n_sig_deg` | FDR < 0.05 且 log2FC > 0.5 的显著 DEG 数量 |
| `mean_log2fc_topN` | Top-N 基因的平均 log2FC |
| `top_genes_str` | Top-20 基因名列表（分号分隔，用于相邻参数组合的 Jaccard 比较） |

#### 热图解读

`param_scan_deg_stability_heatmap.png`：
- 横轴：niche_high quantile（0.70 → 0.85）
- 纵轴：Leiden resolution（0.3 → 0.7）
- 颜色深度：显著 DEG 数量（越深 = 越多显著基因）
- 每格中标注具体数字

**如何使用热图选参数：**

1. 找到颜色最深的区域（DEG 数量最多）
2. 若该区域与相邻格（±0.1 resolution 或 ±0.05 quantile）数值相差不大，说明此区域是稳定的高原区，选取该区域内任一参数组合均合理
3. **若当前参数（resolution=0.5, quantile=0.80）接近热图中颜色最深的区域**：现有参数被证明合理，无需调整
4. **若最优区域落在其他参数**：建议迁移至该参数重跑主分析

**典型结果期望：** Leiden resolution 在 0.4–0.6、quantile 在 0.75–0.80 时通常能得到最多的显著 DEG，因为这个范围产生 6–10 个合理的邻域聚类，能较好区分免疫抑制 niche 与其他区域。

---

## 3. run_chc23_validation.py 所有输出

> 运行命令：`python code/run_chc23_validation.py`（无需参数，使用默认路径）

此脚本完成三个阶段的 CHC23 独立验证，所有输出在 `results/spatial_niche_chc23/` 和 `results/chc23_validation/` 下。

---

### 3.1 CHC23 平行 Niche 分析结果

**目录：** `results/spatial_niche_chc23/`

与 CHC20 主分析（`results/spatial_niche/`）结构完全一致，使用相同参数（k=15, resolution=0.5, quantile=0.80）。

| 文件 | 说明 |
|------|------|
| `spatial_niche_scores.csv` | CHC23 每个 spot 的完整分析数据（与 CHC20 格式相同） |
| `spatial_niche_parameters.csv` | 分析参数记录（含 n_niche_high 等） |
| `immunosuppressive_niche_signature_genes_ranked.csv` | CHC23 Top 特征基因表（含 FDR、pvalue 新列） |
| `immunosuppressive_niche_signature_genes.txt` | CHC23 特征基因列表（纯文本，供 TCGA 分析） |
| `prior_gene_set_auc.csv` | CHC23 先验基因集 AUC 检验结果 |
| `gini_score_genes.csv` | CHC23 Gini Index 评分基因 |
| `sensitivity_analysis.csv` | CHC23 敏感性分析（Jaccard 稳定性） |
| `plots/` | CHC23 所有可视化图（与 CHC20 一一对应） |

#### 如何与 CHC20 对比

将 `results/spatial_niche/plots/` 和 `results/spatial_niche_chc23/plots/` 中同名图片并排比较：
- `spatial_niche_semantic_labels.png`：免疫抑制 niche 的空间分布是否在切片中存在
- `spatial_hepatocyte.png`：肿瘤实质区域分布
- `region_score_boxplots.png`：各区域评分分布

---

### 3.2 跨切片一致性量化结果

**目录：** `results/chc23_validation/`

#### 3.2.1 一致性汇总报告

**文件：** `results/chc23_validation/cross_slice_consistency_report.csv`

| 指标名 | 含义 | 验收标准 |
|--------|------|---------|
| `jaccard` | CHC20 与 CHC23 Top-50 特征基因的 Jaccard 相似度 | **≥ 0.4** |
| `jaccard_pass` | 是否达到验收标准 | `True` / `False` |
| `n_overlap` | 两切片共有特征基因数量 | ≥ 20（对应 Jaccard ≥ 0.4） |
| `spearman_r` | 细胞类型均值比例向量的 Spearman 相关系数 | 越接近 1.0 越好，建议 > 0.7 |
| `spearman_p` | Spearman 检验 p 值 | < 0.05 为显著 |
| `pct_chc20` | CHC20 niche_high spot 占比 | - |
| `pct_chc23` | CHC23 niche_high spot 占比 | - |
| `pct_diff` | 两切片占比差异 | **< 0.05（5%）为一致** |
| `pct_consistent` | 占比一致性是否通过 | `True` / `False` |

#### 如何使用报告

全部三项验收标准（Jaccard ≥ 0.4 + Spearman > 0.7 + 占比差 < 5%）均通过，表明：
> 在 CHC20 上发现的免疫抑制空间生态位在独立的 CHC23 切片中同样存在，两切片发现的特征基因高度重叠，细胞组成模式相似，生态位规模一致，**具有跨患者重现性**。

若部分标准未通过：
- **Jaccard < 0.4**：两切片的特征基因差异较大 → 可能是肿瘤异质性，或 CHC23 细胞反卷积质量问题（与早停修复有关）
- **Spearman < 0.7**：细胞组成模式差异大 → 两患者肿瘤微环境存在实质性差异
- **占比差 > 5%**：生态位规模不同 → 可能与肿瘤分期、免疫浸润程度差异有关

---

#### 3.2.2 特征基因重叠韦恩图

**文件：** `results/chc23_validation/signature_gene_overlap_venn.png`

展示 CHC20 Top-50 和 CHC23 Top-50 特征基因的集合重叠关系：
- 左圆（蓝色）：CHC20 专有基因
- 右圆（橙色）：CHC23 专有基因
- 交叉区域：两切片共有基因（标注具体基因名，最多 15 个）
- 图题标注 Jaccard 系数

#### 3.2.3 细胞类型比例 Spearman 相关图

**文件：** `results/chc23_validation/celltype_spearman_correlation.png`

散点图，每个点代表一种细胞类型：
- 横轴：该细胞类型在 CHC20 所有 spot 中的均值比例
- 纵轴：该细胞类型在 CHC23 所有 spot 中的均值比例
- 对角线：完美相关参考线（y = x）
- 图题：Spearman r 和 p 值

若点集中在对角线附近，说明两切片的细胞组成模式相似（同一类型的肿瘤）。

#### 3.2.4 共同基因表达热图

**文件：** `results/chc23_validation/gene_overlap_heatmap.png`（需要 AnnData 文件可读时才生成）

展示两切片共有特征基因（最多 30 个）在四组 spot（CHC20_High / CHC20_Low / CHC23_High / CHC23_Low）中的 Z-score 标准化均值表达量。

- 每行为一个共有基因，每列为一个 niche 分组
- Z-score 颜色编码（RdBu_r 配色：红色 = 高表达，蓝色 = 低表达）
- 若 CHC20_High 和 CHC23_High 颜色模式一致（均为红色），说明跨切片重现性良好

#### 3.2.5 基因并排列表图

**文件：** `results/chc23_validation/gene_overlap_parallel_list.png`

左右两列分别列出 CHC20 和 CHC23 的 Top-50 特征基因，**共有基因用红色粗体标注**，专有基因用黑色细体。

适合直接放入论文补充材料，直观展示跨切片基因重叠情况。

---

### 3.3 CHC23 DE 验证结果

**文件：** `results/spatial_signature_genes_chc23.txt`

使用 CHC23 自己的 `niche_high` 标签作为分组变量，通过 scanpy Wilcoxon 检验完成差异表达分析后的 Top 特征基因纯文本列表（每行一个基因）。

可直接用于：
1. 与 `results/spatial_signature_genes.txt`（CHC20 主分析结果）对比
2. 作为 TCGA 生存分析的输入（`code/tcga_survival_analysis.R`）
3. 评估两患者的 signature gene 重叠度

---

*文档生成日期：2026-06-25*
