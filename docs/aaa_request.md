# 需求文档

根据以下需求==添加修改代码==

1. 针对CHC23切片regression拟合有问题的情况，后续cell2location的结果并不好，所有的概率基本上都是一个值，箱线图变成了一条线，因此我认为它的拟合结果并不够好，并没有真的收敛，使用你给出的解决方案

   **情况二：确认假收敛（需要调参）**

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


2. CHC23 切片的后续应用，单独写一个代码文件添加一个step来进行验证分析。要求该代码文件中需要有每个函数的注释，并在文件最开头写一个跟其他代码文件类似的介绍。

   **详细方案**

   **阶段一：CHC23 Step 2 重跑（Niche 分析独立验证）**

   直接以 CHC23 反卷积结果为输入，运行与主分析完全相同的参数：

   ```bash
   python code/run_spatial_niche_analysis.py \
       --adata results/adata_vis_post_CHC23.h5ad \
       --out-dir results/spatial_niche_chc23 \
       --signature-out results/spatial_signature_genes_chc23.txt
   ```

   输出与 CHC20 的 `results/spatial_niche/` 结构完全平行，便于逐图对比。

   **阶段二：跨切片一致性量化（需新增分析代码）**

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

   **阶段三：CHC23 Step 3 DE 验证**

   使用 CHC23 自己的 `niche_high` 标签作为分组变量，重跑差异表达分析：

   ```bash
   python code/run_de_analysis.py \
       --adata results/adata_vis_post_CHC23.h5ad \
       --coloc results/spatial_niche_chc23/spatial_niche_scores.csv \
       --coloc-column niche_high \
       --out results/spatial_signature_genes_chc23.txt
   ```

   最终在论文或报告中，将 CHC20 和 CHC23 的 Top-50 signature gene 列表并排展示，突出重叠基因，即构成跨切片重现性的核心证据。

3. Leiden 分辨率与分位数阈值（0.80、0.75）的最优参数组合问题

   **方法二：DEG 稳定性扫描（更严格的客观标准）**

   **核心逻辑**：好的参数组合应使 niche_high vs niche_low 的差异表达结果更显著、更稳定。

   具体操作：

   1. 对 resolution ∈ {0.3, 0.4, 0.5, 0.6, 0.7} 和 niche_high_quantile ∈ {0.70, 0.75, 0.80, 0.85} 构成的所有参数组合（5 × 4 = 20 种）分别运行 niche 分析；

   2. 对每种参数组合，提取 Top-50 signature genes，做 Wilcoxon 检验，记录：
      - FDR < 0.05 的基因数量（显著 DEG 数）
      - Top-50 基因的平均 log2FC

   3. 绘制热图（resolution × quantile），颜色编码 "显著 DEG 数"；

   4. 选择显著 DEG 最多且与相邻参数组合的 Top-50 gene Jaccard 相似度最高（稳健性最好）的参数组合。

   **结论判断标准**：若当前参数（resolution=0.5, quantile=0.80）已经接近热图中的最优区域，则现有参数被证明合理；若最优区域落在其他参数，则迁移到该参数组合。

4. 敏感性分析结果全部相同问题的改进方案

   不比较 spot 数量（永远相同），而是比较**哪些 spot 被选中**是否一致。以主分析（k=15）的 `niche_high` 集合为参考，计算其他参数设置下的 `niche_high` 集合与它的 Jaccard 相似度：

   ```python
   ref_mask = (score_k15 >= score_k15.quantile(0.80))
   
   for each param setting:
       curr_mask = (score_curr >= score_curr.quantile(0.80))
       jaccard = (ref_mask & curr_mask).sum() / (ref_mask | curr_mask).sum()
       # Jaccard 越高，说明两种参数选出的 spot 越重叠，结果越稳健
   ```

   **结果解读**：如果所有参数的 Jaccard 都 > 0.85，说明 niche 的空间位置非常稳健，不依赖邻域参数选择。

5. 筛选出目标免疫抑制基因（FOXP3、TGFB1、FAP 等）

   解决方案：以下三种方案并行

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

6. L-R 通讯图基因只有 5 个且有 2 个全为 0 的修复方案

   **改进 product score 计算方式**

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

   **在缺失基因时寻找同家族替代基因**

   对缺失的基因（如 CCR4），在代码中配置同家族候选基因作为备用：

   ```python
   GENE_FALLBACKS = {
       "CCR4":   ["CCR2", "CCR5"],    # 趋化因子受体家族
       "IL10RA": ["IL10RB"],          # IL-10 受体
   }
   ```

   当主基因缺失时，自动尝试 fallback 基因，并在标签中注明替代。

7. 目前的immunosuppressive niche 空间图与聚类图完全相同，添加基于 niche_high 分位数阈值切割的 spot 空间二值分布图。

   需要新增一张图：**基于 niche_high 分位数阈值切割的 spot 空间二值分布图**，颜色编码为"评分高于 80 百分位（niche_high=True）"与"其余 spot"，从而直观展示"用评分方法选出的区域"与"用聚类方法标注的区域"之间的差异（可能有重叠，也可能有分歧，这个分歧本身就是有信息量的）。

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

8. 将修改过的代码新生成的results结果的详细讲解添加到 ./docs/result_explain.md 中