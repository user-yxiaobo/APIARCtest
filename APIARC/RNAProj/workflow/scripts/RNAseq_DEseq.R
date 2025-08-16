#!/usr/bin/env Rscript

if(!require("yaml")) install.packages("yaml", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("optparse")) install.packages("optparse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
library(DESeq2)
library(yaml)
library(optparse)

.libPaths(c("/usr/local/R/lib64/R/library", .libPaths()))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to the gene count matrix CSV file"),
  make_option(c("-c", "--config"), type = "character", help = "Path to the config JSON string"),
  make_option(c("-o", "--output"), type = "character", help = "Path to the output result CSV file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查是否提供了必要的命令行参数
if (is.null(opt$input) || is.null(opt$config) || is.null(opt$output)) {
  stop("Error: Input, config, and output file paths must be provided. Use --input, --config, and --output options.")
}

cat("[INFO] 输入文件路径: ", opt$input, "\n")
cat("[INFO] 配置文件路径: ", opt$config, "\n")
cat("[INFO] 输出文件路径: ", opt$output, "\n")

# 从指定路径读取基因计数矩阵
cat("[INFO] 正在读取基因计数矩阵...\n")
database_all <- read.table(file = opt$input, sep = ",", header = TRUE)
cat("[INFO] 读取的基因计数矩阵内容如下:\n")
print(head(database_all))

# 从指定路径解析 JSON 配置文件
cat("[INFO] 正在解析配置文件...\n")
config <- yaml.load(opt$config)
cat("[INFO] 配置内容:\n")
print(config)

# 提取表达数据（基因计数矩阵）
cat("[INFO] 正在处理表达矩阵...\n")
expdata <- database_all[, -1]  # 去掉第一列 gene_id
rownames(expdata) <- database_all[, 1]  # 设置行名为基因 ID
cat("[INFO] 表达矩阵的维度: ", dim(expdata), "\n")

# 提取样本名和条件
cat("[INFO] 正在提取样本名称与条件...\n")
sample_names <- names(config)
conditions <- factor(unlist(config), levels = c("T", "P"))
cat("[INFO] 样本名称:\n")
print(sample_names)
cat("[INFO] 样本条件:\n")
print(conditions)

# 检查样本名称与表达矩阵列名是否匹配
if (!all(sample_names %in% colnames(expdata))) {
  cat("[ERROR] 样本名称与表达矩阵列名不匹配！\n")
  cat("表达矩阵的列名:\n")
  print(colnames(expdata))
  cat("配置中的样本名称:\n")
  print(sample_names)
  stop("样本名称与列名不匹配，请检查输入数据和配置文件！")
}

# 构建 colData 数据框
cat("[INFO] 正在构建 colData 数据框...\n")
coldata <- data.frame(conditions = conditions, row.names = sample_names)
cat("[INFO] colData 内容:\n")
print(coldata)
expdata <- expdata[, sample_names]

# 创建 DESeq 数据集
cat("[INFO] 正在移除含 NA 值的行...\n")
expdata <- expdata[complete.cases(expdata), ]
cat("[INFO] 创建 DESeq 数据集...\n")
dds <- DESeqDataSetFromMatrix(countData = expdata, colData = coldata, design = ~ conditions)

# 运行 DESeq2 分析
cat("[INFO] 正在运行 DESeq2 分析...\n")
dds <- DESeq(dds)

# 获取结果并合并标准化表达数据
cat("[INFO] 正在处理分析结果...\n")
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
row.names(resdata) <- resdata[, 1]

names(resdata)[names(resdata) == "Row.names"] <- "SYMBOL"
resdata$gene_id <- sub("\\|[^|]*$", "", resdata$SYMBOL)
resdata$SYMBOL <- sub("^.*\\|", "", resdata$SYMBOL)
resdata <- resdata[!is.na(resdata$padj), ]

cat("[INFO] 分析结果内容如下:\n")
print(head(resdata))

# 保存结果到指定路径
cat("[INFO] 正在保存结果到输出文件...\n")
write.csv(resdata, opt$output, row.names = FALSE)

cat("[INFO] 分析完成，结果已保存到: ", opt$output, "\n")

# Rscript RNAseq_DEseq.R --input /home/yxiaobo/RCProj/RNAseq/database/Chr_enrich/GSE140552/mRNA/gene_count_matrix.csv --config "{'mmChr_enrich_re1': 'T', 'mmChr_enrich_re2': 'T', 'mmNMuMG_veh_re1': 'P', 'mmNMuMG_veh_re2': 'P'}" --output /home/yxiaobo/RCProj/RNAseq/database/Chr_enrich/GSE140552/mRNA/mmNMuMG_veh_vs_Chr_enrich.csv

