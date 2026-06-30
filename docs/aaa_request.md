# 需求文档

仔细阅读我的项目代码，根据以下需求修改代码，删除没有必要的代码，并同步修改readme和每个代码文件最前面的介绍

1. 我目前的项目代码在 Visium spot 级别做 Mann-Whitney U / Wilcoxon rank-sum 时，p 值和 FDR 极小，将这个检验筛选 signature genes 的过程修改为

   ```
   保留条件：
       FDR_block < 0.1
       and AUC_block > 0.60
       and (log2FC > 0.3 or delta_frac > 0.10)
   
   排序：
       signature_rank_score =
           0.35 * minmax(log2FC)
         + 0.35 * minmax(delta_frac)
         + 0.30 * minmax(AUC)
   ```

   在修改保留阈值时判断一下保留的基因个数，如果小于80个的话就放松这个保留阈值。immunosuppressive_niche_signature_genes_ranked.csv 这个表格和同名的txt文件也同步更新。

   将筛选出来的基因过滤掉非特异性基因，例如

   ```
   housekeeping genes
   ribosomal genes
   mitochondrial genes
   hemoglobin genes
   cell-cycle genes
   broad inflammation genes
   tissue damage/stress genes
   ```

   过滤后的所有基因作为后续CHC20验证和大队列生存分析的依据。

2. CHC20 空间转录组数据的验证过程，目前采取的是

   ```
   score_i = mean(log-normalized expression of signature genes)
   ```

   不太合理，结果也不理想。将其修改为将所有保留基因分为三组*（目前没有运行得到的 signature genes 的基因结果，无法作为准确计算，你只需要完整这个结构，保证其能够正常运行，具体三组基因包含哪些我自己修改）*，按照以下规则进行 score 的计算和验证

   三套 score：

   ```
   1. immune_signature_score
   2. stromal_ECM_signature_score
   3. immune_stromal_niche_score
   ```

   具体公式可以这样：

   对每个 spot 或患者，先算 immune 模块：

   ```
   immune_score_i =
       mean(E_ig for g in immune genes)
   ```

   再算 stromal/ECM 模块：

   ```
   stromal_score_i =
       mean(E_ig for g in stromal/ECM genes)
   ```

   然后标准化：

   ```
   z_immune_i = z(immune_score_i)
   z_stromal_i = z(stromal_score_i)
   ```

   最后组合成 niche score：

   ```
   immune_stromal_niche_score_i =
       z_immune_i + z_stromal_i
   ```

   最后这个过程中生成的图片改成：

   ```
   Fig.5A: CHC20 immune module spatial score
   Fig.5B: CHC20 stromal/ECM module spatial score
   Fig.5C: combined immune-stromal niche score
   Fig.5D: HCC4R vs CHC20 combined score distribution
   ```