# 安装Bioconductor和BioStrings包（如果尚未安装）
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

# 加载Biostrings包
library(Biostrings)

# 定义引物序列
primer <- "TCAATATGCTGAAACGCGCGAGAAACCG"
# 定义要在其中搜索的DNA序列（将'DNA.str'替换为您的实际序列）
primer.str <- DNAString(primer)
# 读取基因组数据库文件
sequences <- readDNAStringSet(file.choose(),
                              "fasta")
# 定义一个函数来搜索模式，并在找不到匹配时自动尝试反向互补序列
search_pattern <- function(primer, sequences, max.mismatch, min.mismatch, fixed) {
  matches <- vmatchPattern(primer,
                           sequences,
                           max.mismatch = max.mismatch, 
                           min.mismatch = min.mismatch,
                           fixed = fixed)
  if (length(matches) == 0) {
    rc_primer <- reverseComplement(primer)
    matches <- vmatchPattern( rc_primer, 
                              sequences, 
                              max.mismatch = max.mismatch, 
                              min.mismatch = min.mismatch,
                              fixed = fixed)
  }
  return(matches)
}



# 在DNA序列中搜索模式
results <- search_pattern(primer.str, 
                          sequences, 
                          max.mismatch = 3,
                          min.mismatch = 0,
                          fixed = FALSE)

# 打印结果
print(results)
