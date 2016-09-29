###############################################
###       Compare essential genes
###############################################

library(pROC)
library(GGally)
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

max.old <- read.table("ctx_demo_genetable_collapsed.old.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
max.old <- max.old[-(1:3), c(1,14,13)]    # select ave_gamma, delete neg_control
colnames(max.old) <- c("Gene", "ave_gamma", "p.value")
nrow(max.old)     # num of raw genes: 31977
max.old <- max.old[!grepl("pseudo", max.old$Gene, fixed = TRUE),]
nrow(max.old)     # num of real genes: 15977
max.old$ave_gamma <- as.numeric(max.old$ave_gamma); max.old$p.value <- as.numeric(max.old$p.value)
sum(max.old$Gene %in% sci.K562$Gene)      # num of overlap genes with sci.K562: 15596

CD437 <- read.table("CD437_genetable_collapsed.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
CD437 <- CD437[-c(1:3), c(1,4,3)]   # select gamma, delete neg_control
colnames(CD437) <- c("Gene", "gamma", "p.value")
nrow(CD437)     # num of raw genes: 31954
CD437 <- CD437[!grepl("pseudo", CD437$Gene, fixed = TRUE),]
nrow(CD437)     # num of real genes: 15977
CD437$gamma <- as.numeric(CD437$gamma); CD437$p.value <- as.numeric(CD437$p.value)

### 2. Plot ROC
###############################
# set the sci.K562 p-value < 0.05 to be the true positive, and plot ROC

# exploratory analysis
sci.K562 <- cbind(sci.K562, K562.sig = as.numeric(sci.K562$K562.adjusted.p.value < 0.05))
sum(sci.K562$K562.sig); sum(sci.K562$K562.sig[sci.K562$K562.CS < 0])     # num of hits/positive hits: 2335/1663

max.new <- cbind(max.new, max.new.hit.prob = 1 - max.new$p.value)
sum(max.new$max.new.hit.prob > 0.95); sum(max.new$max.new.hit.prob[max.new$ave_gamma < 0] > 0.95)   # num of possible hits/negtive hits: 2362/1969

max.old <- cbind(max.old, max.old.hit.prob = 1 - max.old$p.value)
sum(max.old$max.old.hit.prob > 0.95); sum(max.old$max.old.hit.prob[max.old$ave_gamma < 0] > 0.95)   # num of possible hits/negtive hits: 2381/1987

CD437 <- cbind(CD437, CD437.hit.prob = 1 - CD437$p.value)
sum(CD437$CD437.hit.prob > 0.95); sum(CD437$CD437.hit.prob[CD437$gamma < 0] > 0.95)   # num of possible hits/negtive hits: 2834/2370

# select negative score
sci.K562 <- sci.K562[sci.K562$K562.CS < 0,]
max.new <- max.new[max.new$ave_gamma < 0,]
max.old <- max.old[max.old$ave_gamma < 0,]
CD437 <- CD437[CD437$gamma < 0,]

# select overlap genes
pred.roc <- merge(sci.K562[,c("Gene", "K562.sig")], max.new[,c("Gene", "max.new.hit.prob")], by = "Gene")
pred.roc <- merge(pred.roc, max.old[,c("Gene", "max.old.hit.prob")], by = "Gene")
pred.roc <- merge(pred.roc, CD437[,c("Gene", "CD437.hit.prob")], by = "Gene")

# plot ROC
png("ROC.05.png", width = 900, height = 700, units = "px")
plot.roc(pred.roc$K562.sig, pred.roc$max.new.hit.prob, col = "red", main = "ROC for different CRISPR screen (using p-value < 0.05 as cut off)")
plot.roc(pred.roc$K562.sig, pred.roc$max.old.hit.prob, col = "blue", add = TRUE)
plot.roc(pred.roc$K562.sig, pred.roc$CD437.hit.prob, col = "green", add = TRUE)
legend("bottomright", legend = c("max.new (AUC:0.729)", "max.old (AUC:0.746)", "CD437 (AUC:0.725)"), col = c("red", "blue", "green"), lty = 1, lwd = 2)
dev.off()


# ### 3. Using more strict cut off
# ######################################
# # set the sci.K562 adjusted p-value < 0.001 to be the true positive, and plot ROC
# sci.K562 <- cbind(sci.K562, K562.sig.adjust = as.numeric(sci.K562$K562.adjusted.p.value < 0.001))
# 
# pred.roc.strict <- merge(sci.K562[,c("Gene", "K562.sig.adjust")], max.new[,c("Gene", "max.new.hit.prob")], by = "Gene")
# pred.roc.strict <- merge(pred.roc.strict, max.old[,c("Gene", "max.old.hit.prob")], by = "Gene")
# pred.roc.strict <- merge(pred.roc.strict, CD437[,c("Gene", "CD437.hit.prob")], by = "Gene")
# 
# plot.roc(pred.roc$K562.sig, pred.roc$max.new.hit.prob, col = "red", main = "ROC for different CRISPR screen (using p-value < 0.001 as cut off)")
# plot.roc(pred.roc$K562.sig, pred.roc$max.old.hit.prob, col = "blue", add = TRUE)
# plot.roc(pred.roc$K562.sig, pred.roc$CD437.hit.prob, col = "green", add = TRUE)
# 
# # no significant improvement in ROC

### 3. Correlation analysis
#########################################

# p.value correlation
p.value <- merge(sci.K562[,c("Gene","K562.adjusted.p.value")], max.new[,c("Gene","p.value")], by = "Gene")
p.value <- merge(p.value, max.old[,c("Gene","p.value")], by = "Gene")
p.value <- merge(p.value, CD437[,c("Gene","p.value")], by = "Gene")
colnames(p.value)[3:5] <- c("Max.new.p.value", "Max.old.p.value", "CD437.p.value")

# ggpairs(p.value[,-1])       # a mess

# # select small p.value to analyze correlation
# p.value <- p.value[p.value$K562.adjusted.p.value < 0.05,]
# ggpairs(p.value[,-1])   # still mess


### 4. Select overlap gene with small p.value
#############################################
p.value <- subset(p.value, K562.adjusted.p.value < 0.001 & Max.new.p.value < 0.001 & Max.old.p.value < 0.001 & CD437.p.value < 0.001)
hit.gene <- p.value$Gene
hit.gene.table <- sci.K562[sci.K562$Gene %in% hit.gene, -4]
hit.gene.table <- merge(hit.gene.table, max.new[,-4], by = "Gene")
hit.gene.table <- merge(hit.gene.table, max.old[,-4], by = "Gene")
hit.gene.table <- merge(hit.gene.table, CD437[,-4], by = "Gene")
colnames(hit.gene.table)[4:9] <- c("max.new.gamma", "max.new.p.value", "max.old.gamma", "max.old.p.value", "CD437.gamma", "CD437.p.value")

write.csv(hit.gene.table, file = "essential.genes.csv", row.names = FALSE)
