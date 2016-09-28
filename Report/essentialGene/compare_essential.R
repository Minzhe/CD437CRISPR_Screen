###############################################
###       Compare essential genes
###############################################

library(pROC)
setwd("~/Project/CD437CRISPR_Screen/Report/essentialGene/")

### 1. read and clean data
################################
sci.K562 <- read.table("Science_essentialCS.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sci.K562 <- sci.K562[,c(1,5,6)]     # select K562
nrow(sci.K562)    # num of genes: 18166

max.new <- read.table("ctx_demo_full_genetable_collapsed.new.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
max.new <- max.new[-(1:3), c(1,14,13)]    # select ave_gamma, delete neg_control
colnames(max.new) <- c("Gene", "ave_gamma", "p.value")
nrow(max.new)     # num of raw genes: 31954
max.new <- max.new[!grepl("pseudo", max.new$Gene, fixed = TRUE),]
nrow(max.new)     # num of real genes: 15977
max.new$ave_gamma <- as.numeric(max.new$ave_gamma); max.new$p.value <- as.numeric(max.new$p.value)
sum(max.new$Gene %in% sci.K562$Gene)      # num of overlap genes with sci.K562: 15596



### 2. Plot ROC
###############################
# set the sci.K562 p-value < 0.05 to be the true positive, and plot ROC

# exploratory analysis
sci.K562 <- cbind(sci.K562, K562.sig = as.numeric(sci.K562$K562.adjusted.p.value < 0.05))
sum(sci.K562$K562.sig); sum(sci.K562$K562.sig[sci.K562$K562.CS < 0])     # num of hits/positive hits: 2335/1663

max.new <- cbind(max.new, max.new.hit.prob = 1 - max.new$p.value)
sum(max.new$max.new.hit.prob > 0.95); sum(max.new$max.new.hit.prob[max.new$ave_gamma < 0] > 0.95)   # num of possible hits/negtive hits: 2362/1969

# select overlap genes
pred.roc <- merge(sci.K562[,c("Gene", "K562.sig")], max.new[,c("Gene", "max.new.hit.prob")], by = "Gene")

# plot ROC
plot.roc(pred.roc$K562.sig, pred.roc$max.new.hit.prob)
