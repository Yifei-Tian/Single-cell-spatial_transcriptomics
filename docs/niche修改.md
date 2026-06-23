# niche修改

**最推荐方案：把当前“手动阈值 niche”改成“文献支持的邻域组成聚类 + 功能评分注释”**

也就是：

1. 对每个 spot/cell 构建空间邻域
   用固定半径或 kNN。比起 `neighbor_radius_multiplier=1.25`，更推荐做成：
   - Visium spot：用空间网格邻接或固定半径；
   - 单细胞空间数据：用 kNN，例如 10-30 个邻居；
   - 同时做 sensitivity analysis：k=10, 20, 30 或 radius=1.0, 1.25, 1.5 倍 spot distance。
2. 计算每个邻域的细胞组成向量
   例如：
   `Hepatocyte, CAF, Endothelial, T cell, Treg, Macrophage, B cell...` 的局部比例。
3. 对这些组成向量做聚类
   用 Leiden / k-means / Gaussian mixture 都可以。这个就是 Schürch 等 “cellular neighborhoods” 的核心思想：不是先定阈值，而是让局部细胞组成自己形成 neighborhood/niche 类型。
4. 再用功能分数给 niche 命名
   例如：
   - `immunosuppressive_score`
   - `Treg_score`
   - `checkpoint_score`
   - `CAF_score`
   - `hypoxia_score`
   - `angiogenesis_score`
5. 最后定义目标 niche
   例如把某个 cluster 注释为：
   “hepatocyte-associated immunosuppressive stromal niche”，依据是：
   - Hepatocyte 邻近度高；
   - stromal/CAF 成分高；
   - immune suppressive module 高；
   - 空间上靠近肿瘤或肝细胞区域。

这样比现在的：

```
hep_high_quantile = 0.75
niche_high_quantile = 0.80
stroma_thr = 0.60
```

更有文献支撑，也更容易在论文/报告里解释。

**可以引用的文献方向**

1. **Schürch et al., Cell, 2020**
   题目：*Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front*
   DOI: [10.1016/j.cell.2020.07.005](https://doi.org/10.1016/j.cell.2020.07.005)
   用途：支持“基于局部细胞组成定义 cellular neighborhoods/niches”。
2. **Palla et al., Nature Methods, 2022, Squidpy**
   题目：*Squidpy: a scalable framework for spatial omics analysis*
   DOI: [10.1038/s41592-021-01358-2](https://doi.org/10.1038/s41592-021-01358-2)
   用途：支持空间邻接图、邻域富集、空间 omics 标准分析流程。你的 Python 项目很适合接 Squidpy。
3. **Keren et al., Cell, 2018**
   题目：*A structured tumor-immune microenvironment in triple negative breast cancer revealed by multiplexed ion beam imaging*
   DOI: [10.1016/j.cell.2018.08.039](https://doi.org/10.1016/j.cell.2018.08.039)
   用途：支持肿瘤-免疫空间结构、局部微环境分型。
4. **CellCharter / NicheCompass / BANKSY 这类方法**
   如果你想更“方法学先进”，可以考虑这些。但对你当前项目来说，我不建议一上来就换成很复杂的模型。你现在的脚本更适合先升级为：
   **Squidpy 邻域图 + composition clustering + module score annotation**。

**你报告里可以这样写**

> 本研究采用“组成为主、功能为辅”的空间 niche 识别框架。首先基于空间邻接关系构建每个 spot/cell 的局部邻域，并计算邻域内主要细胞类型组成；随后对邻域组成向量进行无监督聚类，以识别空间 cellular neighborhoods。最后结合免疫抑制相关基因模块评分、基质细胞比例以及肝细胞邻近度，对候选 niche 进行功能注释。该策略参考了 Schürch et al. 对 cellular neighborhoods 的定义思路，并结合 Squidpy 中空间邻域分析框架。对于邻域半径和分位数阈值，本研究不将其作为固定生物学常数，而是作为数据驱动参数，并通过敏感性分析评估结果稳定性。

**我对你当前参数的具体建议**

- `neighbor_radius_multiplier=1.25`：可以保留，但不要说它来自文献；改成“主分析参数”，并补充 `1.0, 1.25, 1.5` 敏感性分析。
- `hep_high_quantile=0.75`：建议改为注释用，不作为 niche 发现的第一步。
- `niche_high_quantile=0.80`：可以用于选取 high-score niche，但要报告 0.75/0.80/0.85 是否稳定。
- `stroma_thr=0.60`：这个最容易被质疑。建议改成 stromal composition 的连续变量，或者用聚类后再判断哪个 cluster 是 stroma-rich。

一句话：**评价是合理的，但真正该改的不是去找这些阈值的“出处”，而是把阈值规则降级为注释/筛选，把 niche 发现主体改成邻域组成聚类。** 这会更像正式空间转录组分析，也更容易被文献支撑。