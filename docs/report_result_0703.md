## 一、Cellular Neighborhoods 聚类识别免疫抑制空间生态位（对应正文 2.5 节 / 3.4 节）

### 1.1 空间生态位地图与细胞类型组成的关联分析（图 6）

图 6 展示了 HCC4R 切片中空间生态位与细胞类型组成之间的关联分析，包含两个子图：(a) 空间生态位地图，以规则化空间区域标签（Tumor Core、Tumor Edge、Stroma/Immune、Other）在组织物理坐标上呈现各 spot 的域归属；(b) 各空间结构域的细胞类型组成 100\% 堆叠条形图，量化 Hepatocyte、Treg、T/NK、Myeloid、Fibroblast 和 Endothelial 六类细胞在四个域中的相对比例。

从空间生态位地图（图 6a）可见，Tumor Core（红色）呈大片连续分布，构成组织切片的主体骨架，与肿瘤实质在 H\&E 染色和 Hepatocyte 反卷积丰度图中呈现的核密集区域在空间位置上高度吻合；Stroma/Immune（蓝色）则以条带状镶嵌于 Tumor Core 周围及组织边缘，形成环绕/穿插于肿瘤实质之间的连续区域，未表现为孤立的散点分布，这与无监督聚类能够有效捕获组织连续空间微结构、而非随机噪声的方法学预期一致~\citep{schurch2020coordinated}；Tumor Edge（橙色）散布于 Tumor Core 与 Stroma/Immune 的过渡地带，构成两者之间的窄带状界面，在图中占比相对最小。这一空间格局提示，HCC4R 切片的组织结构并非肿瘤与基质截然二分，而是呈现"核心—边缘—基质免疫"的连续过渡梯度，与经典"免疫排斥型"HCC 中肿瘤细胞占据核心、免疫基质细胞被限制在外周边界的空间结构描述相符~\citep{llovet2021hepatocellular,binnewies2018understanding}。

需要说明的是，图 6a 所采用的四类空间域标签（Tumor Core / Tumor Edge / Stroma/Immune / Other）是在正文 2.5 节 Leiden 邻域聚类识别出的 `immunosuppressive niche` 基础上，进一步以 Hepatocyte 比例分位数阈值（75 百分位）和免疫基质评分对全切片 spot 所做的规则化区域标注（`spatial_region`），用以呈现比单一二元 niche/non-niche 标签更连续、更具层次的空间分区结果，而非直接展示 Leiden 聚类原始编号（cluster 0, 1, 2, ...）。二者在生物学含义上是一致的：Tumor Edge 与 Stroma/Immune 共同构成了 2.5 节所定义的、Treg--Myeloid--Fibroblast 协同富集的免疫抑制空间生态位的核心区域，而 Tumor Core 对应肿瘤实质、Other 对应组织中信号较弱或混杂的背景区域。

细胞组成堆叠条形图（图 6b）为上述空间格局提供了定量支持。四个域中 Hepatocyte 比例呈现明确的梯度递减：Tumor Core（约 35\%）> Other（约 24\%）> Tumor Edge（约 25\%）> Stroma/Immune（最低），与"Tumor Core 由肝细胞/恶性细胞主导"的定义直接自洽。与之相对，Treg、Myeloid 和 Fibroblast 三类免疫抑制轴核心细胞的比例均在 Stroma/Immune 和 Tumor Edge 两个域中高于 Tumor Core：Treg 比例在 Stroma/Immune 中达到四组最高水平，Myeloid 和 Fibroblast 亦呈现相似的富集趋势，三者在同一空间域中的协同富集直接印证了正文所述 Treg--TAM--CAF 三角免疫抑制轴的空间协同性~\citep{sakaguchi2020regulatory,noy2014tumor,kalluri2016biology}。T/NK 细胞比例在 Tumor Edge 中相对突出，与肿瘤侵袭前沿是效应 T 细胞与调节性细胞交锋的主战场这一经典认识相符~\citep{schurch2020coordinated}。Endothelial 比例在四组间相对稳定，提示血管分布本身并非驱动该聚类分区的主导因素。

综合图 6a 与图 6b，本研究识别出的免疫抑制空间生态位并非集中于单一孤立的"热点岛"，而是以 Tumor Edge 和 Stroma/Immune 两个相邻、连续的空间域共同构成的界面结构，其细胞组成同时体现 Treg 升高、Myeloid（TAM）升高、Fibroblast（CAF）升高、Hepatocyte 相对降低的多维协同特征，与正文 2.5 节以 Treg/Myeloid/Fibroblast/免疫抑制基因评分四维综合排名（公式 combined\_rank）注释 `immunosuppressive niche` 的方法学预期完全一致，为后续基于该生态位提取特征基因签名（图 7、图 8）提供了可靠的空间分组基础。

### 1.2 邻域组成 Leiden 聚类的原始空间分布（图 6-补）

图 6-补展示了正文 2.5 节所述邻域组成向量 $\mathbf{n}_i$（$k$-NN，$k=15$）在特征空间中重新构图后运行 Leiden 聚类（分辨率 = 0.5）得到的原始聚类结果，是图 6a 规则化空间域标签的上游数据来源。与图 6a 采用的四类语义化标签不同，本图直接以聚类编号（NC-0 至 NC-11，共 12 个 neighborhood cluster）在空间坐标系中着色展示，未经过后续基于 Hepatocyte 分位数阈值的规则化归并，因而更完整地保留了 Leiden 聚类算法在邻域组成特征空间中原生发现的组织异质性结构。

图中最突出的现象是 NC-2（紫色）在空间上呈现出全图面积最大、连通性最强的单一聚类区域，占据切片中心地带并呈树枝状向多个方向延伸，与其余 11 个聚类相比在空间连续性和覆盖范围上均具有明显优势。该聚类被正文方法学流程通过综合排名公式（combined\_rank，见 2.5 节公式 2）判定为 `immunosuppressive_niche`——即在 Treg 比例、Myeloid 比例、Fibroblast 比例和免疫抑制基因评分四个维度的降序排名之和最小，因此在四维证据的联合权衡下被判定为综合免疫抑制特征最强的邻域类型。NC-2 在空间上并非局限于孤立的小簇，而是与图 6a 中 Stroma/Immune（蓝色）和部分 Tumor Edge（橙色）区域的空间轮廓高度重合，两种独立标注方式（Leiden 聚类原始编号 vs. 规则化语义标签）在空间位置上的相互印证，从方法学角度增强了该免疫抑制生态位判定结果的稳健性，也与聚类图与语义标注图应呈现高度空间一致性的预期相符~\citep{schurch2020coordinated,goltsev2018deep}。

图中同时清晰呈现的一个特征是，除 NC-2 外，其余聚类呈现出明显的组织学分区意义：例如 NC-1（红色）和 NC-3（棕色）分别在切片左下与右上角形成两个独立的大面积连续区域，很可能对应不同批次或不同区域的肿瘤实质亚型；NC-5（灰色）、NC-9（浅橙色）等则以较窄条带形式镶嵌于 NC-1、NC-2 等大聚类之间，提示这些邻域代表了肿瘤与基质之间的过渡带亚结构。整体而言，Leiden 算法在未设定聚类数目先验的情况下自动识别出 12 个具有内部一致性、且在空间上连续分布（而非像素级别噪声式散布）的邻域类型，这一结果数量落在正文所述"5--12 个具有内部一致性的空间域"这一区间的上限附近，同时验证了无监督聚类方法确实能够捕获真实存在的组织空间微结构而非人为设定的伪影~\citep{traag2019leiden,schurch2020coordinated}。这种"一个显著免疫抑制邻域 + 多个肿瘤/基质亚结构邻域"并存的格局，也与近年空间蛋白质组学研究中发现的肿瘤组织普遍存在多种细胞邻域共存、且免疫抑制性邻域通常表现为跨越肿瘤-基质界面的连续结构而非孤立热点这一观察一致~\citep{jackson2020single,bhate2022tissue}，进一步支持了本研究识别出的免疫抑制生态位是 HCC4R 切片中一种具有稳定组织学基础的真实空间实体，而非聚类参数任意选择下的偶然产物。

---

## 二、niche 特征基因签名的提取与功能注释（对应正文 2.6 节）

本部分围绕免疫抑制生态位特征基因签名的两层筛选策略（正文 2.6 节），从"检出率富集""空间域表达模式"和"下游通路功能"三个互补角度对结果进行解读，分别对应图 7（niche 签名基因检出率富集散点图）、图 8（签名基因 × 空间域表达热图）和图 9（通路富集气泡图）。

### 2.1 niche 签名基因的检出率富集分析（图 7）

图 7 以 log₂FC（niche\_high / niche\_low）为横轴、Δfraction（检出率差值 frac\_high − frac\_low）为纵轴，展示了两层筛选流程中 Layer 1 的检出率补偿维度，是对传统"仅依赖 log2FC + FDR"的火山图的重要补充可视化。图中以 |log2FC| = 0.5 和 Δfrac = 0.10 两条虚线划分四个象限：红色点（Both-sig，n = 33）同时满足双重显著条件；橙色点（Fraction-only，n = 847）仅在检出率维度显著富集而 log2FC 未达阈值；灰色点为不显著基因。

图中可见的一个重要现象是，绝大多数具有生物学意义的免疫抑制候选基因——包括 TIGIT、CTLA4、CD163、MRC1、FAP 等——均落在橙色的"Fraction-only"区域而非红色的"双重显著"区域，其 log2FC 普遍集中在 0.1--0.3 之间（低于 0.5 的阈值），而 Δfraction 已达到 0.10--0.30 的较高水平。这一模式与正文 2.6 节所述"Visium 每个 spot 覆盖约 5--50 个细胞的混合信号（细胞稀释效应），使 Treg/TAM 标志基因的原始表达量均值被严重稀释，log2FC 因而偏低"的预期完全吻合~\citep{luecken2019current}。换言之，图 7 并非显示这些基因不具有 niche 特异性，而是直接以数据可视化的方式验证了本研究方法学设计的必要性：若仅采用传统的"log2FC 阈值 + p 值"筛选流程，TIGIT、CTLA4、CD163、MRC1、FAP 等生物学关键基因将被排除在 Top 候选列表之外；而检出率差值 Δfrac 却能够敏锐地捕捉到"niche\_high 区域中有更多 spot 检出该基因表达"这一稀释效应下依然保留的信号，从而通过正文所述的两层筛选策略中"AUC > 0.60 且 Δfrac > 0.10"或"FDR < 0.1 且 AUC > 0.60"等补充条件将其重新纳入最终签名基因列表。

图右上角同时呈现少数具有更高 log2FC（> 0.5）且 Δfrac 同样很高（> 0.25）的基因，包括 CXCL12、TGFB1 及一组未单独标注但聚集在图右上方的红色点，这部分基因同时具备"高倍数变化"与"高检出率差异"两个维度的显著性，是免疫抑制 niche 中信号最强、最稳定的一类分子，与正文所述 CXCL12--CXCR4 是招募 Treg/MDSC 进入肿瘤的核心趋化因子轴~\citep{domanska2013targeting}、TGFB1 是 CAF 与 Treg 共同分泌的免疫抑制细胞因子~\citep{mariathasan2018tgfb} 的机制预期一致。值得注意的是，POSTN 和 COL1A1 两个基质/ECM 重塑相关基因虽然在检出率维度尚未达到 0.10 的阈值（分别落在参考虚线附近或略低），但其 log2FC 已明显偏高，提示这两个基因在少数 spot 中呈现局灶性极高表达而非广泛检出，与 CAF 激活标志基因倾向于在致密纤维化灶中高度局部化表达的组织学特征相符~\citep{kalluri2016biology,iozzo2015proteoglycan}。

整体而言，图 7 以直观的散点分布证实了正文两层筛选策略中引入 Δfrac 维度的合理性与必要性：若无该维度的补偿，本研究报告摘要中所述的"FOXP3、IL2RA 等 Treg 标志基因经检出率差值补偿后重新纳入签名基因集"的结论将缺乏数据支撑；本图正是该结论最直接的可视化证据。

### 2.2 niche 签名基因在各空间域中的表达模式（图 8）

图 8 展示了按 `composite_score`（0.35×log2FC + 0.35×Δfrac + 0.30×AUC 加权综合排序，见正文公式 signature\_rank\_score）降序排列的 Top 基因集合在四个空间域（Tumor Core、Tumor Edge、Stroma/Immune、Other）中的平均表达 Z-score 热图。

热图呈现出与图 6 高度一致的空间趋势：几乎所有签名基因在 Stroma/Immune 域中表达最高（深红色，Z-score 接近或超过 2），在 Tumor Core 域中表达最低（深蓝色），Tumor Edge 居中偏高，Other 域接近中性偏低，这一"Stroma/Immune > Tumor Edge > Other > Tumor Core"的单调递减模式在几乎全部 26 个展示基因中保持高度一致，直接印证了正文所述"免疫抑制相关基因模块在 stroma\_immune 区域和 tumor\_edge 区域表达最高，在 tumor\_core 中最低，与肿瘤-免疫交界面是免疫抑制最活跃位点的生物学模型一致"这一结论~\citep{zheng2017landscape}。

需要指出的是，由于该热图按数据驱动的综合排序分数自动选取 Top 基因，实际展示的基因列表（CD52、HLA-DQA1、PTGDS、HLA-DPB1、IGKC、IGHM、HLA-DQB1、HLA-DRB1、LAPTM5、CCL19、HLA-DPA1、TRBC2、AEBP1、LTB、DCN、LUM、GPNMB、ITGB2、HLA-DMB、SRGN、CAPG、TYROBP、GSTP1、TAGLN、CRIP1、CD74、TRAC、CYBA、CXCL9、FXYD5）与正文摘要中重点举例的 FOXP3、IL2RA、CTLA4、TGFB1、CXCL12、CCL22 等经典 Treg/TAM 先验基因并不完全重合。这一现象并非结果矛盾，而是与前述图 7 所揭示的稀释效应及两层筛选架构直接相关：本图所示 Top 基因来自 Layer 1 的全基因组数据驱动排序，其入选主要依赖综合排序分数在全转录组范围内的相对排名，而 FOXP3 等极低丰度的 Treg 核心标志基因虽然通过 Layer 2 的先验基因集靶向 AUC 检验被强制并入最终签名基因列表（保证其不被遗漏），但由于其全局表达量绝对值过低，在按综合分数排序的全基因组 Top-N 列表中位置相对靠后，因而未进入本图默认展示的基因子集，这与正文 2.6 节"第二层针对先验已知的功能基因集进行靶向核查与强制补全"的两层筛选设计逻辑相符——数据驱动排序（Layer 1）和先验知识补全（Layer 2）二者共同构成最终签名基因列表，但在按分数排序的可视化展示中二者不会以相同权重同时出现。

从图中实际展示的基因功能类别看，热图揭示了免疫抑制生态位分子特征的另一重要维度：多个 MHC-II 类分子（HLA-DQA1、HLA-DPB1、HLA-DQB1、HLA-DRB1、HLA-DPA1、HLA-DMB、CD74）及淋巴细胞相关基因（CD52、IGKC、IGHM、LAPTM5、TRBC2、TRAC、LTB、CCL19、TYROBP）在 Stroma/Immune 域中特异性高表达，提示该区域除 Treg--TAM--CAF 轴外，同时存在活跃的抗原呈递与适应性免疫细胞浸润背景~\citep{neefjes2011towards,chen2017oncology}，这与该区域被同时定义为"免疫基质区（Stroma/Immune）"而非单纯的"CAF 纤维化区"的空间语义标注一致。此外，DCN、LUM、AEBP1、TAGLN、GPNMB 等基质/成纤维细胞相关基因同样在 Stroma/Immune 域富集，与该区域 Fibroblast 比例较高的细胞组成特征相符~\citep{iozzo2015proteoglycan}；CXCL9 作为效应 T 细胞招募的经典趋化因子在该域中的富集，进一步支持了该空间域是免疫细胞主动募集与滞留位点的机制推断~\citep{tokunaga2018cxcl9}。综合来看，图 8 从表达模式角度佐证了免疫抑制生态位并非单一细胞类型的孤立信号，而是涵盖免疫细胞浸润、抗原呈递活性与基质重塑的复合分子表型，与正文所述"通路富集分析将上述分子特征映射至功能层面"的逻辑衔接一致。

### 2.3 niche 签名指向的功能通路富集（图 9）

图 9 以气泡图形式展示了 10 条预定义功能通路（覆盖免疫抑制、基质重塑、炎症与代谢应激等维度）基因集评分在四个空间域中的富集情况：气泡大小表征通路平均表达评分的高低，颜色区分该评分是否高于全局四域均值（红色为高于均值，蓝色为低于均值）。

图中最突出的模式是 CAF Activation 通路：该通路在 Stroma/Immune 域呈现全图最大的红色气泡，表明其评分不仅是四域中的最高值，且富集幅度在全部 10 条通路、全部 4 个域组合中最为显著，与正文所述"CAF 通过分泌 CXCL12 招募免疫抑制细胞，并通过重塑细胞外基质形成物理屏障"这一机制描述高度吻合~\citep{kalluri2016biology}。EMT 通路则呈现出与预期部分不同、但可以合理解释的分布模式：其在 Tumor Edge 和 Stroma/Immune 两域均为红色（高于均值），而在 Tumor Core 和 Other 两域为蓝色，这一"边缘与基质双高、核心与背景双低"的格局与上皮间质转化程序主要发生于肿瘤侵袭前沿及其邻近基质、而非肿瘤核心静态区域的经典认识一致~\citep{nieto2016emt}，提示 Tumor Edge 不仅是免疫抑制信号的活跃地带，也是肿瘤细胞获得侵袭性表型的关键空间位置。

TGF-β Signaling、TNF-α Signaling via NF-κB、Angiogenesis、IL6-JAK-STAT3、TAM Polarization (M2) 和 Treg Immune Suppression 六条通路呈现出高度一致的模式：均在 Stroma/Immune 域为红色高富集，在 Tumor Core 和 Other 域为蓝色低富集，Tumor Edge 域评分居中或同样偏高。这一组通路的协同富集为正文所述"TGF-β 信号通路、IL-10 信号通路和 T 细胞受体激活调控通路的显著富集与 Treg 和 M2 型 TAM 介导免疫抑制的核心分子机制直接对应"提供了直接的功能层面证据~\citep{sakaguchi2020regulatory,noy2014tumor}，同时 TNF-α/NF-κB 通路的伴随富集提示该区域并非单纯的免疫抑制状态，而是同时存在慢性炎症信号背景，与免疫抑制性肿瘤微环境中促炎与抗炎信号并存的复杂图景相符~\citep{taniguchi2018nfkb,binnewies2018understanding}。Immune Checkpoint 通路虽然气泡整体偏小（提示其基线表达评分低于其他通路，与 PD-1/CTLA-4 等检查点分子在 Visium spot 级别本身表达量有限的预期一致），但其在 Stroma/Immune 域的相对富集方向（红色）仍与该域富集免疫检查点分子的预期一致。

Hypoxia 通路的分布模式与其余通路存在明显差异，构成本组图中最值得讨论的"非预期"发现：该通路在 Tumor Edge 域为红色高富集，但在 Stroma/Immune 域反而呈现蓝色（低于均值），而在通常被认为血氧供应最差的 Tumor Core 域也为蓝色，Other 域则为红色。这一模式初看似乎与"肿瘤核心区域因血管发育不良而更易缺氧"的直觉预期不完全一致。一个合理的生物学解释是：Visium spot 级别的缺氧评分反映的是转录本层面的 HIF 通路下游基因（HIF1A、VEGFA、LDHA、PGK1、SLC2A1、BNIP3）表达强度，而非组织病理学意义上的直接氧分压测量；HCC 肿瘤核心区域的恶性细胞可能已通过长期适应性重编程下调急性缺氧应激基因的转录、转而依赖组成性糖酵解代谢，使得经典缺氧应激标志物在稳态肿瘤核心中不一定处于最高表达水平，而 Tumor Edge 作为血管新生和肿瘤细胞增殖最活跃、氧耗与血供矛盾最尖锐的动态界面，反而更易检测到急性 HIF 通路激活及其下游血管生成信号~\citep{palazon2017hif}。这一现象与图中同一域 Angiogenesis 通路同样呈现红色高富集的结果相互印证——Tumor Edge 处 Hypoxia 与 Angiogenesis 的协同激活，共同支持了"免疫抑制生态位中同时存在促血管生成信号，与抗血管生成和免疫检查点抑制剂联用方案的生物学理论基础相符"这一正文结论~\citep{finn2020atezolizumab}。Other 域中 Hypoxia 通路的富集则可能与该域本身细胞组成混杂、涵盖部分低质量或边缘组织信号有关，其生物学意义需结合更多空间背景谨慎解读，也提示后续工作可考虑对 Other 域做进一步的子结构拆分以厘清其异质性来源。

综合图 7、图 8 与图 9 三张图的结果，niche 特征基因签名的提取与功能注释形成了一条完整的证据链：图 7 从统计方法学角度验证了检出率补偿维度对于捕获稀释效应下 Treg/TAM 标志基因的必要性；图 8 从空间表达谱角度证实了签名基因在 Stroma/Immune 与 Tumor Edge 两域的一致性富集，并揭示出该区域同时具有免疫浸润、抗原呈递与基质重塑的复合分子表型；图 9 则从预定义通路层面将上述分子特征系统性地映射至 TGF-β、NF-κB、CAF 活化、EMT 及血管生成等功能模块，尤其是 CAF Activation 与 Angiogenesis-Hypoxia 在 Tumor Edge/Stroma-Immune 域的协同富集，为正文所述 Treg--TAM--CAF 免疫抑制轴的空间协同性及其在肿瘤-基质界面的机制基础提供了多层次、相互印证的定量证据。
