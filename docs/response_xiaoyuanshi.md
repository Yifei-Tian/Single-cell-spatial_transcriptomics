# 问题解答文档

> 文档对应 `docs/aaa_request.md` 中提出的 8 个技术问题，结合项目现有代码逐一给出分析与最优方案。

---

## Q1：CHC23 验证切片应如何使用？

### 现状诊断

`run_preprocessing.py` 中 CHC23 已完整走完 Cell2location 反卷积流程（Step 1），输出了：

- `results/adata_vis_post_CHC23.h5ad`
- `results/spot_cell_proportion_CHC23.csv`
- `results/cross_slice_comparison/*.png`（均值比例对比图、Treg KDE 图、箱线图）

但 Step 2（`run_spatial_niche_analysis.py`）和 Step 3（`run_de_analysis.py`）**只处理了 CHC20**，CHC23 的验证并没有真正进入下游分析。

### CHC23 的核心价值

CHC23 的定位是**独立重现性验证（independent replication）**：在另一张来自不同患者的切片上重跑相同的分析流程，证明发现的免疫抑制生态位不是 CHC20 的特异性噪声，而是肝癌中普遍存在的空间结构。

### 详细方案

#### 阶段一：CHC23 Step 2 重跑（Niche 分析独立验证）

直接以 CHC23 反卷积结果为输入，运行与主分析完全相同的参数：

```bash
python code/run_spatial_niche_analysis.py \
    --adata results/adata_vis_post_CHC23.h5ad \
    --out-dir results/spatial_niche_chc23 \
    --signature-out results/spatial_signature_genes_chc23.txt
```

输出与 CHC20 的 `results/spatial_niche/` 结构完全平行，便于逐图对比。

#### 阶段二：跨切片一致性量化（需新增分析代码）

仅靠现有的可视化对比图（均值比例条形图）说服力不足，需要量化指标：

1. **Signature gene 重叠度（Jaccard 相似度）**

   分别读取两张切片的 `immunosuppressive_niche_signature_genes_ranked.csv`，取 Top-50 基因列表，计算 Jaccard 系数：

   ```
   Jaccard = |set_CHC20 ∩ set_CHC23| / |set_CHC20 ∪ set_CHC23|
   ```

   **验收标准：Jaccard ≥ 0.4**（即 50 个基因中有 ≥ 20 个重叠）认为一致性良好。

2. **细胞类型比例 Spearman 相关**

   对每种细胞类型，计算两个切片各 spot 的全局均值比例向量之间的 Spearman 相关系数，评估组成模式的相似程度。

3. **niche_high 占比一致性**

   比较两个切片的 `niche_high` spot 占比（`spatial_niche_parameters.csv` 中 `n_niche_high / n_spots`），如果比例接近（差异 < 5%），说明生态位规模一致。

#### 阶段三：CHC23 Step 3 DE 验证

使用 CHC23 自己的 `niche_high` 标签作为分组变量，重跑差异表达分析：

```bash
python code/run_de_analysis.py \
    --adata results/adata_vis_post_CHC23.h5ad \
    --coloc results/spatial_niche_chc23/spatial_niche_scores.csv \
    --coloc-column niche_high \
    --out results/spatial_signature_genes_chc23.txt
```

最终在论文或报告中，将 CHC20 和 CHC23 的 Top-50 signature gene 列表并排展示，突出重叠基因，即构成跨切片重现性的核心证据。

---

## Q2：CHC23 RegressionModel 只跑了 90 轮提前退出的原因与解决方案

### 根因分析

首先明确一个关键事实：**RegressionModel 训练的输入是 scRNA-seq 数据，不是 Visium spot**。CHC23 spot 数量少于 CHC20 这件事，与 RegressionModel 的训练过程**完全无关**，修改 `batch_size` 不会解决任何问题。

真正的触发机制来自 `preprocessing.py` 中的早停配置：

```python
# preprocessing.py L325-333
early_stopping_patience=30,       # 连续 30 轮无改善则停止
early_stopping_min_delta=1e-4,    # 改善幅度低于此值视为无改善
early_stopping_monitor="elbo_train"
```

CHC23 训练时使用的是 `adata_sc_chc23 = adata_sc[:, shared_genes_chc23]` 这个子集。如果 CHC23 的 Visium 数据质控过滤后保留的基因数（即 `shared_genes_chc23`）少于 CHC20，scRNA 子集的基因维度更小，ELBO 在低维空间中收敛更快，30 轮内没有继续改善，触发早停，这是**正常行为**。

### 如何判断是否真的收敛

查看 `results/regression_training_history_CHC23.png`：

- **真收敛（无需处理）**：曲线在 60 轮后已平滑趋近水平，90 轮时近乎零斜率。
- **假收敛（需要处理）**：曲线在 90 轮时仍有明显下降趋势，或者持续震荡，被过早打断。

### 解决方案

**情况一：确认真收敛（推荐直接接受，无需修改）**

Cell2location 官方文档指出，RegressionModel 的 ELBO 在数百到数千个细胞的 scRNA 数据上通常 100-150 轮内收敛。90 轮正在合理区间内。

==**情况二：确认假收敛（需要调参）**==

在 `run_preprocessing.py` 中针对 CHC23 单独传入更保守的早停参数：

```python
model_chc23 = setup_and_train_regression_model(
    adata_sc_chc23,
    regression_model_cls,
    max_epochs=400,
    early_stopping_patience=50,      # 从 30 增大到 50，给更多观察窗口
    early_stopping_min_delta=5e-5,   # 从 1e-4 降低到 5e-5，对微小改善更敏感
)
```

**为何不应修改 batch_size**：

`batch_size` 控制的是每个 epoch 中处理多少个细胞（scRNA 的 obs），影响的是每轮训练的噪声程度和内存占用，不影响 ELBO 最终收敛到的值域。将 batch_size 从 1024 改小，只会让每轮梯度更嘈杂，反而可能导致训练更不稳定。

---

## Q3：Leiden 分辨率与分位数阈值（0.80、0.75）是否最优？如何寻找最优值？

### 背景说明

代码中涉及两个关键阈值：

- `leiden_resolution=0.5`（当前值）：决定 neighborhood cluster 数量
- `niche_high_quantile=0.80`（当前值）：决定哪些 spot 进入 niche_high 分组（用于 signature 提取）
- `hep_high_quantile=0.75`（当前值）：决定哪些 spot 被标注为 tumor_core（辅助区域注释）

### 方法一：统计学密度分布曲线寻找自然拐点

**用于确定 niche_high_quantile**

绘制 `immunosuppressive_niche_score` 的核密度估计（KDE）曲线，寻找**双峰之间的谷底**作为自然分割点。如果 KDE 曲线呈现明显双峰形态，则谷底对应的分位数即为最自然的阈值（不依赖主观设定的 0.80）。

判断依据：

- 若 score 分布呈**单峰**（接近正态），则 0.80 分位数是合理的，因为只有约 20% 的 spot 属于 niche_high，这与 Visium 切片中免疫抑制微环境通常占少数的生物学预期一致。
- 若 score 分布呈**双峰**，谷底位置对应的分位数（例如 0.72 或 0.85）是更具生物学意义的分割点。

对 `hep_high_quantile` 同样可以绘制 Hepatocyte 比例的 KDE 曲线，找谷底确认 tumor_core 边界。

**用于确定 leiden_resolution**

绘制分辨率扫描图：横轴为 resolution ∈ {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}，纵轴为该分辨率下产生的 neighborhood cluster 数量。当 cluster 数量增长趋于平稳（出现"肘部"）时对应的分辨率即为合理选择。通常 6-10 个 cluster 对于 Visium 切片是生物学合理的范围。

### ==方法二：DEG 稳定性扫描（更严格的客观标准）==

**核心逻辑**：好的参数组合应使 niche_high vs niche_low 的差异表达结果更显著、更稳定。

具体操作：

1. 对 resolution ∈ {0.3, 0.4, 0.5, 0.6, 0.7} 和 niche_high_quantile ∈ {0.70, 0.75, 0.80, 0.85} 构成的所有参数组合（5 × 4 = 20 种）分别运行 niche 分析；

2. 对每种参数组合，提取 Top-50 signature genes，做 Wilcoxon 检验，记录：
   - FDR < 0.05 的基因数量（显著 DEG 数）
   - Top-50 基因的平均 log2FC

3. 绘制热图（resolution × quantile），颜色编码 "显著 DEG 数"；

4. 选择显著 DEG 最多且与相邻参数组合的 Top-50 gene Jaccard 相似度最高（稳健性最好）的参数组合。

**结论判断标准**：若当前参数（resolution=0.5, quantile=0.80）已经接近热图中的最优区域，则现有参数被证明合理；若最优区域落在其他参数，则迁移到该参数组合。

---

## Q4：敏感性分析结果全部相同的原因与改进方案

### 根本原因分析

查看当前敏感性分析代码（`run_spatial_niche_analysis.py` L 603-619），`niche_high` 的定义是：

```python
thr = float(score.quantile(niche_high_quantile))
n_high = int((score >= thr).sum())
```

这里 `niche_high_quantile = 0.80` 是固定值，因此：

- **无论邻域参数 k 或 radius 如何变化，`n_niche_high` 永远等于 `0.20 × n_spots`（约 20% 的 spot）**；
- `score_mean ≈ 0`（Z-score 之和的均值必然接近 0）；
- 结果中所有行的 `n_niche_high = 224` 和 `niche_pct ≈ 20.04%` 完全相同，完全无法区分参数优劣。

这是一个**设计缺陷**：当前敏感性分析只是改变了 score 的计算方式，但最终用固定分位数切割，掩盖了参数差异。

### 改进方案

==**改进一（你的思路：比较 spot 重合度——推荐）**==

不比较 spot 数量（永远相同），而是比较**哪些 spot 被选中**是否一致。以主分析（k=15）的 `niche_high` 集合为参考，计算其他参数设置下的 `niche_high` 集合与它的 Jaccard 相似度：

```python
ref_mask = (score_k15 >= score_k15.quantile(0.80))

for each param setting:
    curr_mask = (score_curr >= score_curr.quantile(0.80))
    jaccard = (ref_mask & curr_mask).sum() / (ref_mask | curr_mask).sum()
    # Jaccard 越高，说明两种参数选出的 spot 越重叠，结果越稳健
```

**结果解读**：如果所有参数的 Jaccard 都 > 0.85，说明 niche 的空间位置非常稳健，不依赖邻域参数选择。

**改进二（分数切割——你的思路之二）**

放弃固定分位数阈值，改为**固定绝对分数阈值**（如 score > 0）进行切割，这样不同参数产生的 score 分布差异就能体现在 `n_niche_high` 的变化上。但这要求先确定一个有意义的 score 零点（Z-score 加总时零点即为"平均水平"），0 分以上即为"高于平均的免疫抑制区域"。

**改进三（DEG 指标）**

对每种参数设置下的 niche_high/niche_low 分组，做 Wilcoxon 检验，记录 FDR < 0.05 的 DEG 数量。DEG 最多的参数组合是"生物学信号最强"的最优邻域参数。

---

## Q5：为什么要选出免疫抑制生态位 cluster？它在后续有什么应用？

### 为什么要选它

免疫抑制生态位（immunosuppressive niche）在肿瘤微环境中代表了一种**功能性协同抑制单元**：

- **Treg 细胞**通过 CTLA4/IL-10/TGF-β 抑制效应 T 细胞（CD8+）的杀伤功能；
- **肿瘤相关巨噬细胞（TAM，Myeloid 细胞）**通过 M2 极化释放 IL-10/TGF-β，同时上调 PD-L1；
- **肿瘤相关成纤维细胞（CAF，Fibroblast 细胞）**构建致密的胶原屏障，物理隔绝效应 T 细胞进入肿瘤核心区域。

三种细胞在空间上的**共定位**是免疫逃逸效率最高的配置。寻找这个 cluster 的目的是精确标注"免疫治疗最难奏效的空间区域"。

### 下游应用（按流程顺序）

1. **Step 3：作为 DE 分析的分组标签**
   `niche_semantic_label == "immunosuppressive_niche"` 的 spot 集合，直接作为 `run_de_analysis.py` 中 Wilcoxon 检验的阳性组，提取该生态位特异性高表达的基因。

2. **Step 4-5：signature gene 投影到 TCGA 做生存分析**
   从这些 spot 中提取的特征基因（如 FOXP3、TGFB1、FAP 等），通过 ssGSEA 投影到 TCGA-LIHC 的 bulk RNA-seq 数据，计算每位患者体内这一"免疫抑制空间程序"的激活程度评分，再做 Cox 回归，证明它与总生存期（OS）显著负相关。这是将空间发现转化为临床预后 biomarker 的关键桥梁。

3. **未来应用（超出当前代码范围）**
   - 预测抗 PD-1/CTLA-4 免疫治疗的疗效（niche 高富集 → 预期疗效差）；
   - 筛选可能逆转这一生态位的药物靶点（如抗 CCL22、抗 CXCL12）；
   - 作为肝癌分型依据（niche-high 型 vs niche-low 型）。

---

## Q6：如何筛选出目标免疫抑制基因（FOXP3、TGFB1、FAP 等）？

### 为什么纯 log2FC 筛不到目标基因

核心原因是**细胞稀释效应**：Visium 每个 spot 包含约 10-30 个细胞，Treg 通常占 spot 的 1-5%，换算下来一个 spot 里平均不到 1 个 Treg 细胞。因此 Treg 特异基因（FOXP3）在 spot 级别的表达量接近 0，log2FC 几乎都是 0，按 log2FC 排序永远排在后面。log2FC 排在前面的往往是在整个 spot 里高表达的肿瘤代谢基因或管家基因。

### 推荐方案（三层递进）

**方案一：Wilcoxon 检验 + 火山图（必做）**

在 `_rank_niche_genes` 函数中引入 p 值计算，将 `scipy.stats.mannwhitneyu` 的 p 值纳入筛选条件：

- 筛选标准：FDR（BH 校正）< 0.05 且 log2FC > 0.5
- 火山图横轴 log2FC，纵轴 -log10(p 值)，设定双阈值虚线，标注目标基因名称
- 这一步能确保筛出的基因有统计显著性支撑，而不仅仅是均值差异

只按 log2FC 排序是当前代码的问题所在，这一步是**修复缺陷**。

**方案二：先验功能基因集检验（针对性发现目标基因）**

将目标基因按生物学功能分为三组，分别对每组做 Mann-Whitney U 检验，报告该基因在 niche_high vs niche_low 中的表达差异是否显著：

```
Treg 标志集: ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"]
TAM 特征集:  ["CD163", "MRC1", "TGFB1", "IL10", "CXCL12"]
CAF 激活集:  ["FAP", "ACTA2", "POSTN", "COL1A1", "CCL22"]
```

对每个基因单独报告：
- 在 niche_high spot 中的均值表达量
- Mann-Whitney U 检验的 p 值（FDR 校正）
- AUC（用于评估区分能力，AUC > 0.6 视为有区分力）

这一方案的价值在于：即使某基因 log2FC 很小（例如 0.2），但如果在 niche_high spot 中的检出率（percent expressed）显著高于 niche_low，p 值仍会显著，而纯 log2FC 会忽略这种差异。

**方案三：Gini Index 特异性评分（补充筛选稀有基因）**

对于 Treg 高度特异的基因（如 FOXP3），其特点是"在少数 spot 中很高，在大多数 spot 中为零"，Gini Index 正好能捕捉这种不均匀性。

对每个目标基因，计算其在 niche_high spot 集合中的 Gini 系数：

```
Gini(x) = 1 - sum((2k - n - 1) * x_k) / (n * sum(x))
# k：排序后的位次
```

Gini 越大（接近 1），说明该基因表达越集中在少数 spot 中，结合 niche_high 标签可以识别高度局灶性表达的免疫抑制基因。

**操作建议**：三种方案并行，取三者都命中的基因作为核心 signature 基因，仅其中一到两种方法命中的基因作为候选基因。

---

## Q7：L-R 通讯图基因只有 5 个且有 2 个全为 0 的原因与修复方案

### 原因分析

LR_PAIRS 中定义了 7 对 L-R 通讯对（`run_spatial_niche_analysis.py` L 100-108）：

```python
LR_PAIRS = [
    ("CXCL12", "CXCR4",   "CXCL12–CXCR4"),
    ("CCL22",  "CCR4",    "CCL22–CCR4"),
    ("TGFB1",  "TGFBR1",  "TGFB1–TGFBR1"),
    ("PDCD1",  "CD274",   "PD-1–PD-L1"),
    ("IL10",   "IL10RA",  "IL10–IL10RA"),
    ("TIGIT",  "NECTIN2", "TIGIT–NECTIN2"),
    ("LAG3",   "HLA-DRA", "LAG3–MHC-II"),
]
```

只显示 5 对说明有 2 对因基因在数据集中不存在而被跳过（代码 L 823-826 的 missing 检查）。2 个全为 0 的原因是：

1. **细胞稀释效应**：与 Q6 同样的问题。CCL22、IL10 等 Treg/TAM 分泌的细胞因子，在 Visium spot 级别几乎检测不到。配体 × 受体的乘积只要任何一方为 0，结果就是 0。

2. **log1p 归一化后的稀疏性**：CP10K + log1p 归一化之后，低表达基因的值虽然不是严格为 0，但均值极小（如 0.0001），乘积几乎为 0，在热图中无法区分。

### 修复方案

**修复一：替换为在 Visium 级别可检测的 L-R 对**

优先保留 CXCL12-CXCR4 和 TGFB1-TGFBR1（这两对通常在 bulk spot 中可检测），补充更多 TAM 相关的高表达通路：

```python
LR_PAIRS_EXPANDED = [
    # 原有可检测的高优先级对
    ("CXCL12", "CXCR4",   "CXCL12–CXCR4"),
    ("TGFB1",  "TGFBR1",  "TGFB1–TGFBR1"),
    ("PDCD1",  "CD274",   "PD-1–PD-L1"),
    # Visium 级别通常可检测的 TAM 通路
    ("SPP1",   "CD44",    "SPP1–CD44"),          # 骨桥蛋白，TAM 分泌
    ("MIF",    "CD74",    "MIF–CD74"),            # 巨噬细胞迁移抑制因子
    ("VEGFA",  "KDR",     "VEGFA–KDR"),           # 血管生成，与 CAF 相关
    ("LGALS9", "HAVCR2",  "Galectin9–TIM-3"),    # 免疫检查点，可检测
    ("CCL2",   "CCR2",    "CCL2–CCR2"),           # 巨噬细胞招募
    # 保留可能有信号的原有对
    ("TIGIT",  "NECTIN2", "TIGIT–NECTIN2"),
    ("LAG3",   "HLA-DRA", "LAG3–MHC-II"),
]
```

==**修复二：改进 product score 计算方式**==

当前计算方式 `product = ligand × receptor` 会因任一方为 0 而清零。改为对配体和受体分别取均值之后再相乘，或者改用加法代替乘法作为通讯强度的代理指标：

```python
# 原有方式（问题在于乘法对稀疏性敏感）
product = expr[ligand].to_numpy() * expr[receptor].to_numpy()

# 改进方式一：对每个 spot 的配体和受体各自做 rank-normalized，再乘
from scipy.stats import rankdata
lig_ranked = rankdata(expr[ligand].to_numpy()) / len(expr)
rec_ranked = rankdata(expr[receptor].to_numpy()) / len(expr)
product = lig_ranked * rec_ranked

# 改进方式二（更常用于 CellChat 等工具的做法）：取配体和受体的几何均值
product = np.sqrt(
    np.maximum(expr[ligand].to_numpy(), 0) *
    np.maximum(expr[receptor].to_numpy(), 0)
)
```

==**修复三：在缺失基因时寻找同家族替代基因**==

对缺失的基因（如 CCR4），在代码中配置同家族候选基因作为备用：

```python
GENE_FALLBACKS = {
    "CCR4":   ["CCR2", "CCR5"],    # 趋化因子受体家族
    "IL10RA": ["IL10RB"],          # IL-10 受体
}
```

当主基因缺失时，自动尝试 fallback 基因，并在标签中注明替代。

---

## Q8：为什么 immunosuppressive niche 空间图与聚类图完全相同？如何添加基于评分的图？

### 原因说明

查看代码 `run_spatial_niche_analysis.py` L 1197-1202：

```python
_spatial_scatter(
    df, "niche_semantic_label",
    plot_dir / "spatial_niche_semantic_labels.png",
    "Immunosuppressive Niche (Leiden Clusters)",
    categorical=True,
)
```

这里 `niche_semantic_label` 列的值直接来源于 `_annotate_neighborhood_clusters` 的输出：其逻辑是把综合排名最低的 Leiden cluster 命名为 `"immunosuppressive_niche"`，其余 cluster 保留为 `"cluster_X"`（L 545-548）。因此这张图本质上就是**Leiden 聚类结果图**，只不过把某一个 cluster 换了个颜色标注出来，与 `spatial_neighborhood_clusters.png` 几乎完全相同，无法体现连续评分的信息。

### 现有评分图

`spatial_immunosuppressive_niche_score.png`（L 1187-1193）已经绘制了连续评分热图，这才是真正"基于分数"的图。但这张图与 niche 发现结果（聚类图）的空间差异无法直观对比。

### 新增方案

==需要新增一张图：**基于 niche_high 分位数阈值切割的 spot 空间二值分布图**==，颜色编码为"评分高于 80 百分位（niche_high=True）"与"其余 spot"，从而直观展示"用评分方法选出的区域"与"用聚类方法标注的区域"之间的差异（可能有重叠，也可能有分歧，这个分歧本身就是有信息量的）。

在 `run_spatial_niche_analysis.py` 的 Step 13（L 1176-1229）中，在 `_spatial_scatter` 的调用序列里**新增以下两个调用**：

```python
# 新增图1：基于评分阈值（niche_high 分位数）选出的 spot 分布图
# 与 spatial_niche_semantic_labels.png 对比，展示聚类与评分方法的异同
_spatial_scatter(
    df, "niche_high",
    plot_dir / "spatial_niche_high_score_spots.png",
    "Niche-High Spots (score > 80th percentile)",
    cmap="RdGy_r",     # 红色=niche_high, 灰色=其他
)

# 新增图2：并排对比图（聚类方法 vs 评分方法）
# 用 matplotlib subplot 将两张图放在同一画面
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
# 左图：Leiden cluster 注释的 niche
# 右图：评分阈值选出的 niche_high spot
# 两图共用同一空间坐标系，颜色说明一致
```

**并排对比图的生物学意义**：

- 如果两张图高度重叠：说明 Leiden 聚类和评分阈值两种方法得到一致结论，结果非常稳健；
- 如果有部分 spot 只在一张图中出现：聚类图中有但评分图中没有的 spot，是"邻域组成像免疫抑制 niche 但自身评分不够高"的 spot（边缘区域）；评分图中有但聚类图中没有的 spot，是"自身评分很高但周围邻域组成不典型"的 spot（孤立的免疫抑制岛）。这两类 spot 的存在是有意义的，分别代表生态位的"社会性"与"自主性"。

---

