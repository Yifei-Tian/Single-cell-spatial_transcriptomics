# ==========================================
# 四、Step 4：TCGA投影（ssGSEA）
# ==========================================
library(TCGAbiolinks)
library(GSVA)

# 假设这里下载并准备好了表达矩阵
# expression matrix
# expr <- ...  # gene x sample

# signature
signature_file <- "../results/spatial_signature_genes.txt"
if (file.exists(signature_file)) {
    geneset <- list("spatial_signature" = scan(signature_file, what=""))

    # 计算 ssGSEA
    # scores <- gsva(expr, geneset, method="ssgsea")

    # 将输出保存
    # write.csv(scores, "../results/tcga_signature_score.csv")
    print("Step 4 完成。")
} else {
    print("未找到空间 signature 基因文件。")
}

# ==========================================
# 五、Step 5：生存分析（核心）
# ==========================================
library(survival)
library(survminer)

# 假设已经将分数与临床数据合并为 df
# df <- merge(score, clinical)

# Cox模型
# cox <- coxph(Surv(time, event) ~ score + age + stage, data=df)
# summary(cox)

# KM曲线
# df$group <- ifelse(df$score > median(df$score), "High", "Low")
# fit <- survfit(Surv(time, event) ~ group, data=df)
# ggsurvplot(fit)

print("Step 5 完成。")

