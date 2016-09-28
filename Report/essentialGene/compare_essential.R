###############################################
###       Compare essential genes
###############################################


### 1. read and clean data
sci.K562 <- read.table("Science_essentialCS.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sci.K562 <- sci.K562[,c(1,5,6)]     # select K562

max.new <- read.table("ctx_demo_full_genetable_collapsed.new.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
max.new <- max.new[-(1:3), c(1,14,13)]    # select ave_gamma, delete neg_control
colnames(max.new) <- c("Gene", "ave_gamma", "p.value")
max.new <- max.new[!grepl("pseudo", max.new$Gene, fixed = TRUE),]
max.new$ave_gamma <- as.numeric(max.new$ave_gamma); max.new$p.value <- as.numeric(max.new$p.value)


