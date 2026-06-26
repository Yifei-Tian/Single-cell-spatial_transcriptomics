# 需求文档

修改代码，并根据修改内容完成代码文件最前面的介绍文字。

1. ### 对于你之前完成的修改

   Layer 1：Wilcoxon + FDR 主签名基因

   **文件路径：** `results/spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv`（格式升级）

   ##### 升级内容

   原文件仅有 `gene / mean_high / mean_low / log2_fc` 四列，新版本新增：

   | 新增列   | 含义                                                     |
   | -------- | -------------------------------------------------------- |
   | `pvalue` | Mann-Whitney U 检验 p 值（单侧，niche_high > niche_low） |
   | `fdr`    | Benjamini-Hochberg 方法校正后的 FDR                      |

   ##### 筛选逻辑变化

   | 筛选维度   | 原版本        | 新版本                                                 |
   | ---------- | ------------- | ------------------------------------------------------ |
   | 统计检验   | 无            | Wilcoxon 秩和检验                                      |
   | 主筛选条件 | 仅 log2FC > 0 | FDR < 0.05 且 log2FC > 0.5                             |
   | 补充条件   | 无            | 若严格条件不足 Top-N，放宽至 FDR < 0.2 且 log2FC > 0.3 |

   ##### 结果解读

   - **FDR < 0.05 的基因**：在 niche_high vs niche_low 中表达差异具有统计学显著性，可安全用于报告
   - **FDR 0.05–0.2 的基因**：有一定差异趋势，作为候选基因，需独立验证
   - **log2FC 含义**：log2FC = 1 意味着 niche_high 的均值表达量是 niche_low 的 2 倍

   ---

   Layer 1 衍生：火山图

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

   ### 存在问题

   多有基因的p值和fdr都很小，而且目标基因也都没有被选择出来，红色的都不是目标基因 

   ### 需要你完成的任务

   解释一下原因，评估找到最佳解决方法，并完成代码修改。以下是供你借鉴的两个策略：

   **策略一：不要只盯着 FC，结合背景表达率（Fraction of spots）** 比起比较总体平均值，比较目标基因在 niche_high 和 niche_low 中的“检出率（表达该基因的 spot 占比）”往往更具有生物学意义。一个典型的微环境标志物可能是：在 niche_high 中有 40% 的 spot 能检测到表达，而在 niche_low 中只有 5% 能检测到。

   **策略二：剥离混杂因素（结合反卷积结果）** 既然你已经做了 Cell2location 得到了各细胞类型的绝对丰度（Abundance），你可以尝试运行一种“基于细胞类型的”差异分析，而不是直接对比原始 spot count。或者将目标区域的 spot 提取出来，单独查看其中基质细胞/免疫细胞 marker 的表达热图，这比粗暴的全图 Wilcoxon 检验更能说明微环境的特性。