# 需求文档

修改并重新整理代码，并根据修改内容完成代码文件最前面的介绍文字。

1. 目前我已经确定HCC4R与CHC20这两个数据集的batch effect非常小，我打算把他们这两个数据集联合分析，修改模式为

   - 不再把其中一个当“验证集”。

   - 用同一个 scRNA 参考只训练一次 `RegressionModel`，得到一套统一的 `cell_state_df`。

   - 用这同一套 reference signature 分别对 `HCC4R` 和 `CHC20` 做 Cell2location 反卷积。

   - 在反卷积之后，把两个切片的 `spot` 级结果合并做联合分析，同时加一个 `sample` 或 `batch` 列标记来源。

   - 下游如果做 cell composition、cluster、niche score、差异比较，就在合并后的表上做；但空间邻域关系最好仍然“按切片内”计算，不要让 HCC4R 的 spot 和 CHC20 的 spot 直接互为物理邻居。

2. 只保留HCC4R和CHC20的数据集联合分析的相关代码，其他不需要的代码删除掉，可选的step3也删除掉
3. 为我的分析添加这两张图，图在 ./docs/plots 中
4. 模仿论文 Single-cell landscape of the ecosystem in earlyrelapse hepatocellular carcinoma 修改生成图片的配色，论文pdf位于 ./docs/Single-cell landscape of the ecosystem in early- relapse hepatocellular carcinoma.pdf
5. 生成的 niche_fraction_scatter.png 和 niche_signature_volcano.png 里面对于基因名字的标注，只标注与超过阈值的和低于阈值但在阈值附近的，其他不显著的重要基因就不要标注了。
6. 根据修改后的内容修改更新 result_explain.md ，使其符合目前的代码版本。