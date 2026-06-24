#!/usr/bin/env Rscript

# ==============================================================================
# 脚本总览：tcga_survival_analysis.R
# ==============================================================================
#
# 【整体任务】
#   本脚本是整个分析流程的第 4-5 步，负责将 Python 端（run_spatial_niche_analysis.py）
#   从空间转录组数据中提炼出的"免疫抑制空间生态位基因签名"，投影到 TCGA-LIHC
#   （肝细胞癌大队列）的 bulk RNA-seq 数据中，计算每位患者的免疫抑制程度评分，
#   再结合临床随访数据进行预后分析。
#
# 【与 Python 分析的衔接关系】
#   Python 端：
#     scRNA-seq + Visium 空间转录组
#       → Cell2location 反卷积（每个空间 spot 的细胞类型丰度）
#       → 构建免疫抑制生态位评分（immunosuppressive_niche_score）
#       → 识别高免疫抑制 spot，比较其与低免疫抑制 spot 的基因差异
#       → 基于 CHC20 主分析输出 Top-80 上调基因 → results/spatial_signature_genes.txt  ← 关键接口
#       → CHC23 作为空间验证集，不直接作为 TCGA 投影输入
#   R 端（本脚本）：
#     读取 spatial_signature_genes.txt
#       → Step 4: ssGSEA 对 TCGA-LIHC ~370 例患者打分
#       → Step 5: KM 曲线 + 多变量 Cox 回归验证预后意义
#
# 【生物学意义】
#   单细胞/空间转录组样本量极小（通常 2-5 例），发现了"肿瘤实质-Treg-髓系细胞"
#   共定位的免疫抑制生态位。但结论是否具有普遍性和临床意义，需要在大队列验证。
#   TCGA-LIHC 提供了 ~370 例肝癌患者的 bulk RNA-seq + 长期随访数据，
#   是验证预后意义的金标准队列。
#   若 Cox 回归显示 HR > 1 且 P < 0.05，说明"免疫抑制生态位特征"在排除年龄、
#   分期等混杂因素后，仍是独立不良预后因子，具有临床转化价值。
#
# 【输入文件】
#   results/spatial_signature_genes.txt  — Python 端（run_spatial_niche_analysis.py）输出的
#                                          免疫抑制生态位特征基因列表（每行一个基因名）
#   （TCGA-LIHC 数据由脚本通过 TCGAbiolinks 包自动从 GDC 数据库下载，无需本地准备）
#     · TCGA-LIHC RNA-seq count 矩阵（~370 例患者）
#     · TCGA-LIHC 临床随访数据（含 OS 时间、生存状态、年龄、分期等）
#
# 【主要输出文件】
#   results/tcga_signature_score.csv      — 每位患者的 ssGSEA 评分
#   results/tcga_signature_survival.csv   — 评分 + 临床信息合并表
#   results/cox_results.txt               — 多变量 Cox 回归摘要
#   results/km_plot.png                   — Kaplan-Meier 生存曲线图
#   results/cox_forest_plot.png           — Cox 回归森林图
# ==============================================================================

# ------------------------------------------------------------------------------
# 加载所需 R 包
# suppressPackageStartupMessages() 用于抑制包加载时的提示信息，保持输出整洁
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(TCGAbiolinks))  # 从 GDC 数据库下载 TCGA 数据
suppressPackageStartupMessages(library(GSVA))          # 提供 ssGSEA 基因集富集打分方法
suppressPackageStartupMessages(library(survival))      # Cox 比例风险模型、KM 曲线核心计算
suppressPackageStartupMessages(library(survminer))     # 生存曲线可视化（ggsurvplot、ggforest）
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required. Please install it with BiocManager::install('SummarizedExperiment').")
}

# ------------------------------------------------------------------------------
# 工具函数 1：get_script_path()
# 【任务】找到本脚本文件的绝对路径，作为解析其他相对路径的基准。
# 【实现】
#   R 脚本被调用的方式有两种：
#     ① 命令行 Rscript 运行：通过 commandArgs() 中的 "--file=..." 参数拿到路径
#     ② RStudio 交互式运行：通过 sys.frames()[[1]]$ofile 拿到路径
#   若两种方式都拿不到，则退而使用当前工作目录 getwd()。
# ------------------------------------------------------------------------------
get_script_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)  # 获取所有命令行参数（含 R 内部参数）
    file_arg <- "--file="
    # 在所有参数中找到以 "--file=" 开头的项，去掉前缀后就是脚本路径
    script_path <- sub(file_arg, "", args[grep(file_arg, args)])
    if (length(script_path) > 0) {
        return(normalizePath(script_path))  # normalizePath() 将路径转为标准绝对路径
    }
    # 若是在 RStudio 中 source() 执行，尝试从调用帧中获取脚本路径
    if (!is.null(sys.frames()[[1]]$ofile)) {
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
    # 兜底方案：返回当前工作目录
    return(normalizePath(getwd()))
}

# ------------------------------------------------------------------------------
# 工具函数 2：parse_args()
# 【任务】解析命令行参数，支持 "--key value" 和 "--key=value" 两种格式，
#         未提供的参数使用 defaults 中的默认值。
# 【实现】
#   1. 若命令行中包含 "--help"，打印使用说明后直接退出。
#   2. 初始化结果 out 为默认值列表。
#   3. 遍历所有参数：
#      - 若参数包含 "="（如 --project=TCGA-LIHC），直接拆分 key 和 value。
#      - 否则（如 --project TCGA-LIHC），将下一个参数作为 value，索引 +2 跳过。
# ------------------------------------------------------------------------------
parse_args <- function(args, defaults) {
    # 若用户传入 --help，打印所有可用参数的说明并退出
    if ("--help" %in% args) {
        cat("Usage: Rscript tcga_survival_analysis.R [--project TCGA-LIHC] [--workflow 'HTSeq - FPKM']\n")
        cat("       [--assay 'HTSeq - FPKM'] [--signature results/spatial_signature_genes.txt]\n")
        cat("       [--out-dir results] [--score-out tcga_signature_score.csv]\n")
        cat("       [--merged-out tcga_signature_survival.csv] [--cox-out cox_results.txt]\n")
        cat("       [--km-out km_plot.png] [--forest-out cox_forest_plot.png] [--seed 1234]\n")
        cat("       [--gdc-retries 3] [--gdc-wait 5] [--gdc-method api] [--files-per-chunk 20] [--min-genes 5]\n")
        cat("       [--ssl-verify 1] [--query-rds path/to/query.rds]\n")
        quit(save = "no", status = 0)
    }
    out <- defaults  # 以默认值为基础，用命令行参数覆盖
    i <- 1
    while (i <= length(args)) {
        item <- args[i]
        if (grepl("^--", item)) {  # 只处理以 "--" 开头的参数
            key <- sub("^--", "", item)  # 去掉 "--" 前缀，得到 key（可能含 "="）
            if (grepl("=", key, fixed = TRUE)) {
                # 格式：--key=value，用 "=" 分割，取左侧为 key，右侧为 value
                kv <- strsplit(key, "=", fixed = TRUE)[[1]]
                out[[kv[1]]] <- kv[2]
                i <- i + 1
                next
            }
            if (i + 1 <= length(args)) {
                # 格式：--key value，下一个参数就是 value，索引跳 2 位
                out[[key]] <- args[i + 1]
                i <- i + 2
                next
            }
        }
        i <- i + 1
    }
    return(out)
}

# ------------------------------------------------------------------------------
# 工具函数 3：resolve_path()
# 【任务】将路径统一转为绝对路径。
# 【实现】
#   判断路径是否以盘符（Windows：C:\）或 "/" 开头：
#     - 是 → 已经是绝对路径，直接标准化返回。
#     - 否 → 相对路径，拼接到 root 目录下再标准化返回。
# ------------------------------------------------------------------------------
resolve_path <- function(path, root) {
    # 判断是否为绝对路径（Windows 盘符 或 Unix "/"）
    if (grepl("^([A-Za-z]:[\\\\/]|/)", path)) {
        return(normalizePath(path, winslash = "\\", mustWork = FALSE))
    }
    # 相对路径：拼接到项目根目录 root 之下
    return(normalizePath(file.path(root, path), winslash = "\\", mustWork = FALSE))
}

# ------------------------------------------------------------------------------
# 工具函数 4：apply_ssl_settings()
# 【任务】当用户指定 --ssl-verify 0 时，关闭 HTTPS 证书验证。
# 【使用场景】在某些网络环境（如企业内网、VPN）下，GDC 服务器的 SSL 证书可能
#             无法通过验证，导致下载失败。关闭验证可以绕过此问题，但有安全风险，
#             仅在信任的网络环境中使用。
# 【实现】
#   通过 httr 包的 set_config() 将 ssl_verifypeer 和 ssl_verifyhost 设为 0，
#   同时设置系统环境变量禁用 curl 的 SSL 吊销检查。
# ------------------------------------------------------------------------------
apply_ssl_settings <- function(ssl_verify) {
    if (!is.na(ssl_verify) && ssl_verify == 0) {
        if (requireNamespace("httr", quietly = TRUE)) {
            # 关闭 httr（R 的 HTTP 客户端）对 SSL 证书的验证
            httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
        }
        # 关闭 curl 级别的 SSL 吊销检查
        Sys.setenv(R_CURL_SSL_REVOKE_BEST_EFFORT = "true")
        message("SSL verification disabled for GDC requests (use only if you trust the network).")
    }
}

# ------------------------------------------------------------------------------
# 工具函数 5：classify_gdc_error()
# 【任务】判断网络错误的类型，以便给出更有针对性的错误提示。
# 【实现】
#   读取错误信息文本，用关键词匹配判断是 SSL 错误、GDC 服务器错误还是其他错误。
# ------------------------------------------------------------------------------
classify_gdc_error <- function(err) {
    msg <- conditionMessage(err)
    if (grepl("SSL|ssl|secure|certificate|connect error", msg)) {
        return("SSL")
    }
    if (grepl("server down|status", msg)) {
        return("GDC server")
    }
    return("Other")
}

extract_valid_workflows <- function(err) {
    msg <- conditionMessage(err)
    lines <- trimws(unlist(strsplit(msg, "\n", fixed = TRUE)))
    candidates <- lines[grepl("^=>", lines)]
    if (!length(candidates)) {
        return(character(0))
    }
    trimws(sub("^=>\\s*", "", candidates))
}

build_gdc_query <- function(project, workflow) {
    GDCquery(
        project       = project,
        data.category = "Transcriptome Profiling",
        data.type     = "Gene Expression Quantification",
        workflow.type = workflow
    )
}

make_gdc_path_component <- function(x) {
    gsub("[^A-Za-z0-9.-]+", "_", x)
}

ensure_gdc_download_dirs <- function(query, directory) {
    results <- getResults(query)
    if (!nrow(results)) {
        return(invisible(FALSE))
    }

    id_col <- intersect(c("file_id", "id"), colnames(results))
    if (!length(id_col)) {
        message("Could not infer GDC file directories from query results; continuing without pre-creating UUID folders.")
        return(invisible(FALSE))
    }

    project <- if ("project" %in% colnames(results)) unique(results$project) else character(0)
    if (!length(project) || is.na(project[1])) {
        project <- query$project
    }
    project <- project[1]

    data_category <- if ("data_category" %in% colnames(results)) unique(results$data_category)[1] else "Transcriptome Profiling"
    data_type <- if ("data_type" %in% colnames(results)) unique(results$data_type)[1] else "Gene Expression Quantification"
    data_category <- make_gdc_path_component(data_category)
    data_type <- make_gdc_path_component(data_type)

    target_dirs <- file.path(
        directory,
        project,
        data_category,
        data_type,
        results[[id_col[1]]]
    )
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    invisible(vapply(target_dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

resolve_assay_config <- function(se, requested_assay, workflow) {
    assays <- SummarizedExperiment::assayNames(se)
    if (!length(assays)) {
        stop("No assay found in GDCprepare result.")
    }

    selected_assay <- requested_assay
    if (identical(requested_assay, "auto")) {
        candidates <- if (identical(workflow, "STAR - Counts")) {
            c("tpm_unstrand", "fpkm_unstrand", "fpkm_uq_unstrand", "unstranded")
        } else {
            assays
        }
        selected_assay <- candidates[candidates %in% assays]
        selected_assay <- if (length(selected_assay)) selected_assay[1] else assays[1]
    }

    if (!(selected_assay %in% assays)) {
        message(sprintf("Assay '%s' not found; using '%s'", selected_assay, assays[1]))
        selected_assay <- assays[1]
    }

    counts_like <- grepl("count|unstranded|stranded", selected_assay, ignore.case = TRUE) &&
        !grepl("fpkm|tpm", selected_assay, ignore.case = TRUE)

    list(
        assay_name = selected_assay,
        kcdf = if (counts_like) "Poisson" else "Gaussian",
        log_transform = !counts_like
    )
}

# ------------------------------------------------------------------------------
# 工具函数 6：gdc_retry()
# 【任务】带自动重试机制地执行 GDC 数据库操作，应对网络不稳定问题。
# 【实现】
#   采用"指数退避"重试策略：
#     - 第 1 次失败后等待 wait × 1 秒再重试
#     - 第 2 次失败后等待 wait × 2 秒再重试
#     - 以此类推，直到达到最大重试次数
#   若全部重试耗尽仍失败，调用 classify_gdc_error() 判断错误类型，
#   给出具体的排错建议（如 SSL 问题建议 --ssl-verify 0）后抛出错误终止脚本。
# 【参数说明】
#   fn      : 待执行的函数（如 GDCquery / GDCdownload / GDCprepare）
#   retries : 最大重试次数（默认 3）
#   wait    : 基础等待时间（秒，实际等待 = wait × 尝试次数）
#   label   : 操作名称，用于打印提示信息
# ------------------------------------------------------------------------------
run_ssgsea <- function(expr, geneset, kcdf, verbose = FALSE) {
    expr <- as.matrix(expr)

    if (exists("ssgseaParam", envir = asNamespace("GSVA"), mode = "function")) {
        param_fun <- get("ssgseaParam", envir = asNamespace("GSVA"))
        param_args <- list(exprData = expr, geneSets = geneset)
        param_formals <- names(formals(param_fun))
        if ("alpha" %in% param_formals) {
            param_args$alpha <- 0.25
        }
        if ("normalize" %in% param_formals) {
            param_args$normalize <- TRUE
        }
        param <- do.call(param_fun, param_args)
        return(tryCatch(
            GSVA::gsva(param, verbose = verbose),
            error = function(e) GSVA::gsva(param)
        ))
    }

    GSVA::gsva(
        expr,
        geneset,
        method       = "ssgsea",
        kcdf         = kcdf,
        abs.ranking  = TRUE,
        verbose      = verbose
    )
}

make_finite_hr <- function(x, lower = 1e-3, upper = 1e3, missing = 1) {
    x <- as.numeric(x)
    x[is.na(x)] <- missing
    x[is.infinite(x) & x > 0] <- upper
    x[is.infinite(x) & x < 0] <- lower
    pmin(pmax(x, lower), upper)
}

save_cox_forest_plot <- function(cox_fit, data, path) {
    tryCatch(
        {
            forest_plot <- survminer::ggforest(cox_fit, data = data)
            ggplot2::ggsave(path, plot = forest_plot, width = 6, height = 5, dpi = 150)
            message(sprintf("Forest plot saved: %s", path))
        },
        error = function(e) {
            message(sprintf("ggforest failed: %s", conditionMessage(e)))
            message("Creating a capped fallback Cox forest plot instead.")

            fit_summary <- summary(cox_fit)
            ci <- as.data.frame(fit_summary$conf.int)
            coef_table <- as.data.frame(fit_summary$coefficients)
            if (!nrow(ci)) {
                stop("No Cox coefficients available for forest plot.", call. = FALSE)
            }

            p_col <- grep("^Pr\\(", colnames(coef_table), value = TRUE)
            p_value <- if (length(p_col)) coef_table[[p_col[1]]] else NA_real_
            plot_df <- data.frame(
                term = rownames(ci),
                hr = make_finite_hr(ci[["exp(coef)"]]),
                lower = make_finite_hr(ci[["lower .95"]]),
                upper = make_finite_hr(ci[["upper .95"]]),
                p_value = p_value,
                stringsAsFactors = FALSE
            )
            plot_df$term <- factor(plot_df$term, levels = rev(plot_df$term))

            fallback_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(y = term, x = hr)) +
                ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
                ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = term), linewidth = 0.6) +
                ggplot2::geom_point(size = 2.2) +
                ggplot2::scale_x_log10(
                    limits = c(1e-3, 1e3),
                    breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)
                ) +
                ggplot2::labs(
                    x = "Hazard ratio (95% CI, capped to 1e-3 - 1e3)",
                    y = NULL
                ) +
                ggplot2::theme_bw(base_size = 11) +
                ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

            ggplot2::ggsave(path, plot = fallback_plot, width = 6, height = 5, dpi = 150)
            utils::write.csv(plot_df, sub("\\.png$", "_table.csv", path), row.names = FALSE)
            message(sprintf("Fallback forest plot saved: %s", path))
        }
    )
}

gdc_retry <- function(fn, retries, wait, label) {
    for (attempt in seq_len(retries)) {
        # tryCatch 捕获错误：若 fn() 出错，result 存储错误对象而非终止脚本
        result <- tryCatch(fn(), error = function(e) e)
        if (!inherits(result, "error")) {
            return(result)  # 执行成功，直接返回结果
        }
        if (attempt < retries) {
            # 未达最大重试次数，打印进度后等待一段时间再重试
            message(sprintf("%s failed (attempt %d/%d): %s", label, attempt, retries, conditionMessage(result)))
            Sys.sleep(wait * attempt)  # 指数退避：等待时间随失败次数递增
            next
        }
        # 已达最大重试次数，分析错误类型并给出排错提示
        err_type <- classify_gdc_error(result)
        hint <- ""
        if (err_type == "SSL") {
            hint <- "Hint: check HTTPS/SSL settings or try --ssl-verify 0."
        } else if (err_type == "GDC server") {
            hint <- "Hint: GDC may be down; retry later."
        }
        stop(sprintf("%s failed after %d attempt(s): %s. %s", label, retries, conditionMessage(result), hint), call. = FALSE)
    }
    stop(sprintf("%s failed after %d attempt(s).", label, retries), call. = FALSE)
}

# ==============================================================================
# 项目根目录与参数配置
# ==============================================================================

# 获取本脚本的路径，用于定位项目根目录
script_path <- get_script_path()

# 尝试使用硬编码的 Windows 开发路径；若不存在（如在其他机器运行），
# 则自动切换为脚本所在目录的上一级目录作为项目根目录
default_root <- "E:\\Research!!\\Codes\\Bioinformation_train\\Single-cell&spatial_transcriptomics"
project_root <- if (dir.exists(default_root)) {
    normalizePath(default_root, winslash = "\\", mustWork = FALSE)
} else {
    # dirname(script_path) 获取脚本所在目录，".." 向上一级即项目根目录
    normalizePath(file.path(dirname(script_path), ".."), winslash = "\\", mustWork = FALSE)
}

# ------------------------------------------------------------------------------
# 默认参数配置
# 所有参数均可通过命令行覆盖，例如：
#   Rscript tcga_survival_analysis.R --project TCGA-LIHC --seed 42
# ------------------------------------------------------------------------------
defaults <- list(
    project    = "TCGA-LIHC",                              # TCGA 肿瘤队列名称（肝细胞癌）
    workflow   = "STAR - Counts",                           # GDC 数据处理流程类型（FPKM 表达量）
    assay      = "auto",                           # 从 SummarizedExperiment 中读取的 assay 名称
    signature  = file.path("results", "spatial_signature_genes.txt"),  # Python 端基于 CHC20 主分析生成的签名基因列表
    out_dir    = "results",                                 # 所有输出文件的目录
    score_out  = "tcga_signature_score.csv",               # 输出：每位患者的 ssGSEA 评分
    merged_out = "tcga_signature_survival.csv",            # 输出：评分 + 临床信息合并表
    cox_out    = "cox_results.txt",                        # 输出：Cox 回归结果文本
    km_out     = "km_plot.png",                            # 输出：KM 生存曲线图
    forest_out = "cox_forest_plot.png",                    # 输出：Cox 森林图
    seed       = "1234",                                   # 随机数种子（保证结果可重复）
    min_genes  = "5",                                      # 签名基因最少数量（< 5 则报错）
    gdc_retries = "3",                                     # GDC 网络请求最大重试次数
    gdc_wait    = "5",                                     # 重试基础等待时间（秒）
    gdc_method  = "api",                                   # GDCdownload 下载方法：api 或 client
    files_per_chunk = "20",                                # 每个下载分块包含的文件数；越小越稳但越慢
    gdc_dir     = "E:\\GDCdata",                           # Short path avoids Windows path length limits
    ssl_verify  = "1",                                     # SSL 验证开关（1=开启，0=关闭）
    query_rds   = ""                                       # 可选：预先缓存的 GDCquery 结果路径
)

# 解析命令行参数（用户提供的参数覆盖 defaults 中的默认值）
args <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

# 将字符串参数转换为整数，确保后续数值计算正确
args$seed        <- as.integer(args$seed)
args$min_genes   <- as.integer(args$min_genes)
args$gdc_retries <- as.integer(args$gdc_retries)
args$gdc_wait    <- as.integer(args$gdc_wait)
args$files_per_chunk <- as.integer(args$files_per_chunk)
args$ssl_verify  <- as.integer(args$ssl_verify)
if (is.na(args$files_per_chunk) || args$files_per_chunk < 1) {
    stop("--files-per-chunk must be a positive integer.")
}
if (!(args$gdc_method %in% c("api", "client"))) {
    stop("--gdc-method must be either 'api' or 'client'.")
}

# ------------------------------------------------------------------------------
# 路径解析与输出目录创建
# ------------------------------------------------------------------------------

# 将签名文件路径和输出目录解析为绝对路径
signature_path <- resolve_path(args$signature, project_root)
out_dir        <- resolve_path(args$out_dir, project_root)
gdc_data_dir   <- resolve_path(args$gdc_dir, project_root)

# 创建输出目录（recursive=TRUE 支持多级目录，showWarnings=FALSE 避免目录已存在时报警告）
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gdc_data_dir, recursive = TRUE, showWarnings = FALSE)

# 构建各输出文件的完整路径（out_dir 目录 + 文件名）
score_out  <- file.path(out_dir, args$score_out)   # ssGSEA 评分表
merged_out <- file.path(out_dir, args$merged_out)  # 合并后的生存分析数据表
cox_out    <- file.path(out_dir, args$cox_out)     # Cox 回归结果文本
km_out     <- file.path(out_dir, args$km_out)      # KM 曲线图片
forest_out <- file.path(out_dir, args$forest_out)  # 森林图图片

# ------------------------------------------------------------------------------
# 读取签名基因列表（来自 Python 端的输出）
# 这是 Python 分析与 R 分析之间唯一的数据接口
# ------------------------------------------------------------------------------

# 检查签名文件是否存在，不存在则提前报错（避免后续流程白跑）
if (!file.exists(signature_path)) {
    stop(sprintf("Signature file not found: %s", signature_path))
}

# scan() 逐行读取文本文件中的基因名（每行一个基因）
signature_genes <- scan(signature_path, what = "", quiet = TRUE)
signature_genes <- unique(trimws(signature_genes))  # 去重 + 去除首尾空白字符
signature_genes <- signature_genes[nzchar(signature_genes)]  # 去掉空行（nzchar 判断字符串非空）

# 签名基因太少时 ssGSEA 无统计意义，提前终止
if (length(signature_genes) < args$min_genes) {
    stop("Signature gene list too small; please provide at least 5 genes.")
}

# 根据参数决定是否关闭 SSL 验证（网络环境问题时使用）
apply_ssl_settings(args$ssl_verify)

# ==============================================================================
# Step 4: TCGA 数据下载与 ssGSEA 投影打分
# ==============================================================================
# 【任务】
#   从 GDC（Genomic Data Commons）数据库下载 TCGA-LIHC 队列的基因表达数据，
#   对每位患者计算"免疫抑制空间生态位签名"的 ssGSEA 富集分数。
# 【ssGSEA 的生物学意义】
#   ssGSEA 分数高 → 该患者肿瘤中，与空间免疫抑制生态位相关的基因整体表达水平高
#   → 推测该患者肿瘤微环境（TME）中存在类似的免疫抑制空间结构
#   这是"从空间精细解析到大队列 bulk 验证"的标准转化策略。
# ==============================================================================

query <- NULL  # 初始化查询对象为空

# 若用户提供了预缓存的 query RDS 文件（上次运行时保存的），直接加载，跳过网络请求
if (nzchar(args$query_rds)) {
    query_rds_path <- resolve_path(args$query_rds, project_root)
    if (!file.exists(query_rds_path)) {
        stop(sprintf("Query RDS not found: %s", query_rds_path))
    }
    query <- readRDS(query_rds_path)
    message(sprintf("Loaded cached query: %s", query_rds_path))
}

# 若没有缓存，向 GDC 发起查询请求，获取 TCGA-LIHC 的 HTSeq-FPKM 转录组数据信息
if (is.null(query)) {
    message(sprintf("Downloading TCGA expression for %s...", args$project))
    query <- tryCatch(
        gdc_retry(
            function() build_gdc_query(args$project, args$workflow),
            retries = args$gdc_retries,
            wait    = args$gdc_wait,
            label   = "GDCquery"
        ),
        error = function(e) {
            valid_workflows <- extract_valid_workflows(e)
            if (length(valid_workflows) && !(args$workflow %in% valid_workflows)) {
                fallback_workflow <- valid_workflows[1]
                message(sprintf(
                    "Workflow '%s' is unavailable; retrying with '%s'.",
                    args$workflow,
                    fallback_workflow
                ))
                args$workflow <<- fallback_workflow
                return(gdc_retry(
                    function() build_gdc_query(args$project, args$workflow),
                    retries = args$gdc_retries,
                    wait    = args$gdc_wait,
                    label   = sprintf("GDCquery[%s]", args$workflow)
                ))
            }
            stop(e)
        }
    )
}

# 记录下载是否成功的标志，用于最后的状态提示
download_ok <- TRUE
ensure_gdc_download_dirs(query, gdc_data_dir)

# 将查询结果对应的实际数据文件下载到本地
# 若下载失败（如网络中断），不直接终止——后续 GDCprepare() 会尝试使用本地缓存数据
tryCatch(
    {
        gdc_retry(
            function() GDCdownload(
                query,
                directory = gdc_data_dir,
                method = args$gdc_method,
                files.per.chunk = args$files_per_chunk
            ),
            retries = args$gdc_retries,
            wait    = args$gdc_wait,
            label   = "GDCdownload"
        )
    },
    error = function(e) {
        download_ok <<- FALSE  # <<- 修改上层环境的变量（类似全局赋值）
        message(conditionMessage(e))
        message("Download failed; GDCprepare requires all query files, so the analysis will stop before preparation.")
    }
)

if (!download_ok) {
    stop(
        paste(
            "GDCdownload did not complete all files.",
            "Rerun with a smaller chunk size, for example --files-per-chunk 5,",
            "or try --gdc-method client if the GDC Data Transfer Tool is installed."
        ),
        call. = FALSE
    )
}


# 将下载的多个文件整合为一个 SummarizedExperiment 对象（se）
# SummarizedExperiment 是 Bioconductor 的标准数据格式，包含：
#   assay：基因表达矩阵（行=基因，列=样本）
#   rowData：基因注释信息（包含 gene_name）
#   colData：样本/患者临床信息
se <- gdc_retry(
    function() GDCprepare(query, directory = gdc_data_dir),
    retries = args$gdc_retries,
    wait    = args$gdc_wait,
    label   = "GDCprepare"
)

# 检查指定的 assay 名称是否存在，若不存在则自动选第一个可用的 assay
assay_config <- resolve_assay_config(se, args$assay, args$workflow)
assay_name <- assay_config$assay_name
message(sprintf("Using workflow '%s' with assay '%s' (kcdf=%s).", args$workflow, assay_name, assay_config$kcdf))

# ------------------------------------------------------------------------------
# 提取表达矩阵，并将 Ensembl ID 映射到 Gene Symbol
# 【为什么要做这个映射？】
#   TCGA 数据的行名（基因ID）是 Ensembl ID（如 ENSG000...），
#   而签名基因列表使用的是 Gene Symbol（如 FOXP3、CTLA4）。
#   需要将 Ensembl ID → Gene Symbol，才能找到签名基因的表达值。
# ------------------------------------------------------------------------------

# 提取基因表达矩阵（行=基因Ensembl ID，列=样本 TCGA barcode）
expr <- SummarizedExperiment::assay(se, assay_name)

# 从基因注释信息中获取 Gene Symbol（基因通用名称）
# 不同版本的 TCGA 数据列名可能不同，依次尝试几种常见的列名
gene_symbol <- SummarizedExperiment::rowData(se)$gene_name
if (is.null(gene_symbol)) {
    gene_symbol <- SummarizedExperiment::rowData(se)$external_gene_name
}
if (is.null(gene_symbol)) {
    # 若注释表中没有 Symbol 信息，直接用行名（可能是 Ensembl ID）
    gene_symbol <- rownames(expr)
}
gene_symbol <- as.character(gene_symbol)

# 处理缺失的 Symbol：对于 NA 或空字符串的基因，用其原始 Ensembl ID 代替
missing_symbol <- is.na(gene_symbol) | gene_symbol == ""
if (any(missing_symbol)) {
    gene_symbol[missing_symbol] <- rownames(expr)[missing_symbol]
}

# 过滤掉无法获得有效 Symbol 的基因（极少数情况）
valid <- !is.na(gene_symbol) & gene_symbol != ""
expr        <- expr[valid, , drop = FALSE]
gene_symbol <- gene_symbol[valid]

# ------------------------------------------------------------------------------
# 处理重复基因名：将映射到同一 Gene Symbol 的多个 Ensembl ID 的表达值取均值
# 【为什么有重复？】
#   一个 Gene Symbol 可能对应多个 Ensembl ID（如不同转录本），
#   需要合并为一个代表性值，通常取均值。
# ------------------------------------------------------------------------------

gene_table <- table(gene_symbol)  # 统计每个 Symbol 对应多少个 Ensembl ID
# rowsum() 按 gene_symbol 分组对表达矩阵求行和（聚合同名基因）
expr_sum  <- rowsum(expr, group = gene_symbol, reorder = FALSE)
# sweep() 将每行除以该 Symbol 对应的 Ensembl ID 数量，得到均值
expr_mean <- sweep(expr_sum, 1, as.numeric(gene_table[rownames(expr_sum)]), "/")

# 对表达量做 log2(FPKM + 1) 变换：
# +1 是为了避免 log(0)；log2 变换使数据近似正态分布，适合 ssGSEA 的高斯核
expr_for_gsva <- if (assay_config$log_transform) log2(expr_mean + 1) else expr_mean

# ------------------------------------------------------------------------------
# 使用 ssGSEA（Single Sample Gene Set Enrichment Analysis）对每个样本打分
# 【ssGSEA 原理】
#   对每个患者样本，计算签名基因集的富集分数：
#   即签名基因整体表达水平相对于全基因组中所有基因的排序位置。
#   分数越高，表示签名基因集在该样本中整体表达越强。
# 【参数说明】
#   method = "ssgsea"   : 使用单样本 GSEA 方法（每个样本独立计算，不依赖其他样本）
#   kcdf = "Gaussian"   : 核密度估计使用高斯核，适合连续型数据（FPKM），
#                         若是整数型 counts 则应使用 "Poisson"
#   abs.ranking = TRUE  : 使用表达量的绝对排名（而非相对排名）
# ------------------------------------------------------------------------------

geneset <- list(spatial_signature = signature_genes)  # 将签名基因包装为基因集列表
set.seed(args$seed)  # 设置随机数种子确保结果可重复
scores <- run_ssgsea(expr_for_gsva, geneset, assay_config$kcdf, verbose = FALSE)
# scores 是一个 1 行 × N列 的矩阵（1个基因集 × N个样本），提取第一行即全部样本的评分

# ------------------------------------------------------------------------------
# 整理评分结果：从样本级别聚合到患者级别
# 【为什么需要聚合？】
#   TCGA 中一位患者可能有多个样本（如原发肿瘤 + 转移灶），
#   TCGA barcode 前 12 位是患者 ID（如 TCGA-BC-A10Q），
#   取同一患者所有样本的均值，确保每位患者只有一个评分。
# ------------------------------------------------------------------------------

sample_ids  <- colnames(expr_for_gsva)                 # 完整 TCGA 样本条形码（长度 ≥ 15）
patient_ids <- substr(sample_ids, 1, 12)               # 截取前 12 位作为患者 ID

# 构建评分数据框：样本ID、患者ID、ssGSEA 分数
score_df <- data.frame(
    sample_id  = sample_ids,
    patient_id = patient_ids,
    score      = as.numeric(scores[1, ]),              # scores 第一行即 spatial_signature 的分数
    stringsAsFactors = FALSE
)

# 按患者 ID 分组，对多个样本取均值 → 每位患者一个代表性评分
score_patient <- aggregate(score ~ patient_id, data = score_df, FUN = mean)

# 保存评分结果到 CSV 文件
write.csv(score_patient, score_out, row.names = FALSE)
message(sprintf("ssGSEA scores saved: %s", score_out))

# ==============================================================================
# Step 5: 生存分析
# 【任务】
#   将 ssGSEA 评分与 TCGA 患者的临床随访数据合并，
#   通过多变量 Cox 回归和 Kaplan-Meier 曲线验证评分的预后意义。
# 【生物学假说】
#   免疫抑制生态位评分越高 → 肿瘤微环境免疫抑制越重 → 患者预后越差
#   若 Cox 回归 HR > 1 且 P < 0.05，说明该假说在大队列中得到验证。
# ==============================================================================

message("Fetching TCGA clinical data...")

# 从 GDC 获取 TCGA-LIHC 的临床数据（包含患者年龄、分期、生死状态、随访时间等）
clinical <- gdc_retry(
    function() GDCquery_clinic(project = args$project, type = "clinical"),
    retries = args$gdc_retries,
    wait    = args$gdc_wait,
    label   = "GDCquery_clinic"
)

# 兼容不同版本 TCGA 数据中患者 ID 列名的差异
# 新版用 "submitter_id"，旧版用 "bcr_patient_barcode"
if (!("submitter_id" %in% colnames(clinical)) && ("bcr_patient_barcode" %in% colnames(clinical))) {
    clinical$submitter_id <- clinical$bcr_patient_barcode
}
clinical$patient_id <- clinical$submitter_id  # 统一使用 patient_id 作为主键

# ------------------------------------------------------------------------------
# 构建生存分析所需的两个核心变量
# 【生存分析基础概念】
#   time（生存时间）: 患者从诊断/入组到"事件发生"或"失访"的天数
#   event（事件状态）: 1 = 已发生终点事件（死亡）; 0 = 右删失（失访或仍存活）
#
# 【右删失的含义】
#   部分患者随访结束时仍存活（或失联），我们只知道"至少活了这么久"，
#   这种数据称为"右删失"（right-censored），生存分析专门处理这类不完整数据。
# ------------------------------------------------------------------------------

# 生存时间：有死亡记录用死亡时间，否则用最后随访时间（右删失）
clinical$time <- ifelse(
    !is.na(clinical$days_to_death),
    clinical$days_to_death,          # 已死亡：用死亡天数
    clinical$days_to_last_follow_up  # 仍存活/失联：用最后随访天数
)

# 事件标志：有死亡时间 = 1（事件发生），无死亡时间 = 0（删失）
clinical$event <- ifelse(!is.na(clinical$days_to_death), 1, 0)

# 提取年龄（天→年），作为 Cox 回归的协变量（控制年龄混杂效应）
if ("age_at_diagnosis" %in% colnames(clinical)) {
    clinical$age <- clinical$age_at_diagnosis / 365.25  # 将天数转换为年龄（岁）
}

# ------------------------------------------------------------------------------
# 提取肿瘤分期信息，作为 Cox 回归的协变量（控制分期混杂效应）
# TCGA 数据中分期列名可能有差异，依次检查常见列名
# ------------------------------------------------------------------------------

stage_column <- NA
if ("ajcc_pathologic_stage" %in% colnames(clinical)) {
    stage_column <- "ajcc_pathologic_stage"   # 首选：AJCC 病理分期
} else if ("tumor_stage" %in% colnames(clinical)) {
    stage_column <- "tumor_stage"              # 备选：肿瘤分期
}
if (!is.na(stage_column)) {
    clinical$stage <- as.character(clinical[[stage_column]])
    # 将 "Stage I"、"Stage II" 等规范化为 "I"、"II"，去掉 "Stage " 前缀
    clinical$stage <- gsub("^Stage ", "", clinical$stage, ignore.case = TRUE)
}

# ------------------------------------------------------------------------------
# 合并评分数据与临床数据
# 以 patient_id 为主键进行内连接（inner join）：只保留两边都有数据的患者
# ------------------------------------------------------------------------------

df <- merge(score_patient, clinical, by = "patient_id")

# 过滤掉生存时间无效的样本（time 为 NA 或 ≤ 0 的样本无法用于生存模型）
df <- df[!is.na(df$time) & df$time > 0, , drop = FALSE]

# 若有效样本数太少，无法拟合统计模型，提前报错
if (nrow(df) < 5) {
    stop("Not enough samples with survival data to fit models.")
}

# ------------------------------------------------------------------------------
# 构建多变量 Cox 比例风险回归模型
# 【Cox 回归的目的】
#   评估 ssGSEA 评分是否是"独立"的预后因子：
#   即在控制年龄、分期等已知临床因素后，评分是否仍与生存显著相关。
#   若是，说明空间免疫抑制生态位的预后价值不仅仅来自年龄大或分期晚这些混杂因素。
# 【协变量动态添加策略】
#   必含：score（ssGSEA 免疫抑制评分）
#   若数据中有有效年龄信息 → 加入 age 作为协变量
#   若数据中有有效分期信息 → 加入 stage 作为协变量
# ------------------------------------------------------------------------------

# 从必需的 score 开始，根据数据可用性动态添加其他协变量
covariates <- c("score")
if ("age" %in% colnames(df) && any(!is.na(df$age))) {
    covariates <- c(covariates, "age")   # 添加年龄为协变量
}
if ("stage" %in% colnames(df) && any(!is.na(df$stage))) {
    covariates <- c(covariates, "stage")
    df$stage <- as.factor(df$stage)     # 分期转为因子型（Cox 模型需要）
}

# 动态构建 Cox 回归公式，例如：Surv(time, event) ~ score + age + stage
# Surv(time, event) 是生存分析的响应变量，同时编码了时间和事件状态
cox_formula <- as.formula(paste("Surv(time, event) ~", paste(covariates, collapse = " + ")))

# 拟合 Cox 比例风险模型
# coxph() 估计每个协变量的风险比（HR）及其置信区间和 P 值
cox_fit <- coxph(cox_formula, data = df)

# 将 Cox 回归完整摘要（HR、P值、一致性指数 C-index 等）输出到文本文件
sink(cox_out)           # 将后续所有 print() 输出重定向到文件
print(summary(cox_fit)) # summary() 包含 HR、95% CI、Wald检验、Log-rank 检验等
sink()                  # 恢复输出到控制台
message(sprintf("Cox results saved: %s", cox_out))

# ------------------------------------------------------------------------------
# Kaplan-Meier（KM）生存曲线
# 【KM 曲线的目的】
#   直观展示高评分组 vs 低评分组的生存差异，是最易读懂的生存分析可视化方式。
# 【分组策略】
#   以 ssGSEA 评分的中位数为阈值，将患者分为"High"（高免疫抑制）和"Low"两组。
# 【图形含义】
#   - 曲线越快下降，该组患者死亡越快
#   - risk.table = TRUE：在图下方显示各时间点的"at-risk"人数（仍在随访中的人数）
#   - pval = TRUE：在图上显示 log-rank 检验的 P 值（两组生存差异是否显著）
# ------------------------------------------------------------------------------

# 按中位数将患者分为高评分组（High）和低评分组（Low）
df$group <- ifelse(df$score > median(df$score, na.rm = TRUE), "High", "Low")

# 拟合 KM 生存曲线（分组比较）
km_fit <- survfit(Surv(time, event) ~ group, data = df)

# 用 ggsurvplot() 绘制带 at-risk 表和 P 值的 KM 曲线
km_plot <- ggsurvplot(km_fit, data = df, risk.table = TRUE, pval = TRUE)
ggsave(km_out, plot = km_plot$plot, width = 6, height = 5, dpi = 150)
message(sprintf("KM plot saved: %s", km_out))

# ------------------------------------------------------------------------------
# Cox 回归森林图（Forest Plot）
# 【森林图的目的】
#   将多变量 Cox 回归结果可视化：
#   - 每行代表一个协变量（score、age、stage 等）
#   - 水平线段表示该变量的 HR（风险比）及 95% 置信区间
#   - 若置信区间不跨越 HR=1（虚线），说明该变量对预后有显著影响
#   - HR > 1：该变量越高，死亡风险越大（不良预后因子）
#   - HR < 1：该变量越高，死亡风险越小（保护因子）
# ------------------------------------------------------------------------------

save_cox_forest_plot(cox_fit, df, forest_out)

# 保存合并后的完整数据表（评分 + 临床信息 + 分组），供后续分析或审阅
write.csv(df, merged_out, row.names = FALSE)
message(sprintf("Merged survival table saved: %s", merged_out))

# 若数据下载步骤失败但 GDCprepare 使用缓存数据成功，给出友好提示
if (!download_ok) {
    message("GDCprepare succeeded using cached data.")
}

message("Step 4-5 completed.")
