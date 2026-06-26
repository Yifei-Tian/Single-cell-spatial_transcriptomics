# 需求文档

修改代码

1. 目前你给我做了这样的改正
   run_spatial_niche_analysis.py 新增结果

   2.1 敏感性分析（改进版）

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

   | 列名             | 含义                                | 典型值         |
   | ---------------- | ----------------------------------- | -------------- |
   | `param_mode`     | 参数类型（`knn` 或 `radius`）       | `knn`          |
   | `param_label`    | 参数描述                            | `kNN k=10`     |
   | `is_reference`   | 是否为参考参数（k=15）              | `True`/`False` |
   | `n_niche_high`   | niche_high spot 数量                | ~500           |
   | `niche_pct`      | niche_high spot 占比（%)            | ~20.0          |
   | `jaccard_vs_ref` | 与参考集合（k=15）的 Jaccard 相似度 | 0.72–1.00      |
   | `score_std`      | niche 评分的标准差                  | 任意正数       |

   #### 结果解读

   | Jaccard 值 | 含义                                             |
   | ---------- | ------------------------------------------------ |
   | > 0.85     | 参数选出的空间区域与主分析高度一致，**结果稳健** |
   | 0.70–0.85  | 有一定差异，可能存在参数依赖性，建议关注         |
   | < 0.70     | 参数变化对结果影响较大，需重新审视参数选择       |

   这个改正内容有如下几个问题

   1. **损失了连续分布的量级信息 (Magnitude Loss)** Jaccard 是一个无权重的集合指标。它只在乎某个 spot 的评分是否过了 80% 的及格线，完全忽略了评分的绝对大小。
   2. **固定边缘概率导致的数学降维** 由于通过固定分位数强制使得 $|A| = |B| = 224$，Jaccard 公式的分母被死死限制住了。
   3. **边界极度敏感 (Threshold Instability)** 处于 80% 阈值边缘的 spot 会引入巨大的随机噪声。如果一个 spot 在参考组排第 224 名，在实验组排第 225 名（由于极其微小的分数扰动），它就会被踢出集合 $B$，从而导致 Jaccard 下降。这种下降并非因为空间模式发生实质性改变，仅仅是因为人为设定的硬截断（Hard Thresholding）。

   我有这样几个解决方案，你评估一下是否合理，使用最优的方案修改代码的这个部分

   方法一：**Spearman 等秩相关系数 (Spearman's Rank Correlation)**：直接计算不同参数下所有 spot Niche 评分排名的相关性。这保留了所有 spot 的信息，不受 80% 截断值的影响。

   方法二：**连续型 Jaccard (Continuous/Weighted Jaccard)**：如果不截断，而是将每个 spot 的 Niche 评分归一化到 $[0, 1]$ 之间，直接计算连续向量之间的重叠度。

2. 根据目前每个代码文件内容重新修改每个代码文件的开头介绍和总结