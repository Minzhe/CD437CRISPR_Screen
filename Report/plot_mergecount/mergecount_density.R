### 1. Read data ###
library(ggplot2)
setwd("E:/project/CrisprScreen/report/")

## 1. max old data
count.max.old <- read.table("ctx_demo_rawcountstable.old.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
count.max.old <- count.max.old[-(1:3),]
count.max.old[,2:7] <- apply(count.max.old[,2:7], 2, as.numeric)
row.names(count.max.old) <- NULL
count.max.old <- cbind(count.max.old, T0.ave = apply(count.max.old[c("T0", "T0.1")], 1, mean))
count.max.old <- cbind(count.max.old, untreated.ave = apply(count.max.old[c("untreated", "untreated.1")], 1, mean))
# enrichment
count.max.old <- cbind(count.max.old, untre_T0 = count.max.old$untreated.ave - count.max.old$T0.ave)
count.max.old <- cbind(count.max.old, log2enrich.gamma = log2((count.max.old$untreated.ave/sum(count.max.old$untreated.ave))/(count.max.old$T0.ave/sum(count.max.old$T0.ave))))

## 2. max new data
count.max.new <- read.table("ctx_demo_full_rawcountstable.new.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
count.max.new <- count.max.new[-(1:3),]
count.max.new[,2:7] <- apply(count.max.new[,2:7], 2, as.numeric)
row.names(count.max.new) <- NULL
count.max.new <- cbind(count.max.new, T0.ave = apply(count.max.new[c("T0", "T0.1")], 1, mean))
count.max.new <- cbind(count.max.new, untreated.ave = apply(count.max.new[c("untreated", "untreated.1")], 1, mean))
# enrichment
count.max.new <- cbind(count.max.new, untre_T0 = count.max.new$untreated.ave - count.max.new$T0.ave)
count.max.new <- cbind(count.max.new, log2enrich.gamma = log2((count.max.new$untreated.ave/sum(count.max.new$untreated.ave))/(count.max.new$T0.ave/sum(count.max.new$T0.ave))))

## 3. CD437 old data
count.CD437.old <- read.table("CD437_rawcountstable.old.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
count.CD437.old <- count.CD437.old[-(1:3),]
count.CD437.old[,2:4] <- apply(count.CD437.old[,2:4], 2, as.numeric)
row.names(count.CD437.old) <- NULL
# enrichment
count.CD437.old <- cbind(count.CD437.old, untre_T0 = count.CD437.old$untreated - count.CD437.old$T0)
count.CD437.old <- cbind(count.CD437.old, log2enrich.gamma = log2((count.CD437.old$untreated/sum(count.CD437.old$untreated))/(count.CD437.old$T0/sum(count.CD437.old$T0))))




### 2. Plot distribution

# 1. T0 distribution
plot(density(count.max.old$T0.ave), col = "red", lwd = 2, main = "T0 count distribution", xlim = c(-50, 700), ylim = c(0, 0.0065))
lines(density(count.max.new$T0.ave), col = "blue", lwd = 2)
lines(density(count.CD437.old$T0), col = "green", lwd = 2)
legend(x = "topright", legend = c("max.old", "max.new", "CD437"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

# 2. untreated distribution
plot(density(count.max.old$untreated.ave), col = "red", lwd = 2, main = "Untreated count distribution", xlim = c(-50, 700), ylim = c(0, 0.0065))
lines(density(count.max.new$untreated.ave), col = "blue", lwd = 2)
lines(density(count.CD437.old$untreated), col = "green", lwd = 2)
legend(x = "topright", legend = c("max.old", "max.new", "CD437"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

# 3. Untreated - T0 density
untre_T0.old <- data.frame(untre_T0 = count.max.old$untre_T0[abs(count.max.old$untre_T0) > 50], group = "max.old")
untre_T0.new <- data.frame(untre_T0 = count.max.new$untre_T0[abs(count.max.new$untre_T0) > 50], group = "max.new")
untre_T0.CD437 <- data.frame(untre_T0 = count.CD437.old$untre_T0[abs(count.CD437.old$untre_T0) > 50], group = "CD437")
untre_T0 <- rbind(untre_T0.old, untre_T0.new, untre_T0.CD437)
rm(untre_T0.old, untre_T0.new, untre_T0.CD437)
ggplot(untre_T0, aes(untre_T0, fill = group, colour = group)) +
      geom_density(alpha = 0.4, lwd = 0.8, adjust = 0.5) +
      xlim(-400, 200) + ggtitle("Untreadted - T0 (>50)")

# compared unnormalized and normailzed density plot
# ggplot(untre_T0, aes(untre_T0, fill = group, colour = group)) +
#   geom_histogram(breaks = seq(-400, 200, 5), alpha = 0.5, position = "identity", lwd = 0.2) +
#   ggtitle("Untreadted - T0 (>50) unnormalized")
# 
# ggplot(untre_T0, aes(untre_T0, fill = group, colour = group)) +
#   geom_histogram(aes(y = ..density..), breaks = seq(-400, 200, 5), alpha = 0.5, position = "identity", lwd = 0.2) +
#   ggtitle("Untreadted - T0 (>50) Normalized")

### 3. Log2 enrichment
log2enrich <- merge(count.max.old[c("X", "log2enrich.gamma")], count.max.new[c("X", "log2enrich.gamma")], by = "X")
log2enrich <- merge(log2enrich, count.CD437.old[c("X", "log2enrich.gamma")], by = "X")
colnames(log2enrich) <- c("sgID", "max.old", "max.new", "CD437")
log2enrich.melt <- melt(log2enrich, id = "sgID", variable.name = "group")
ggplot(log2enrich.melt[abs(log2enrich.melt$value) > 1,], aes(value, fill = group, colour = group)) +
      geom_density(alpha = 0.4, lwd = 0.8, adjust = 0.5) +
      ggtitle("Log2 enrichment") + xlim(-8, 3)

### 4. Correlation
log2enrich.filter <- log2enrich[!(is.na(log2enrich$max.new) | log2enrich$max.new %in% c(Inf, -Inf)),]
log2enrich.filter <- log2enrich.filter[!(abs(log2enrich.filter$max.old)<0.5 | abs(log2enrich.filter$max.new)<0.5 | abs(log2enrich.filter$CD437)<0.5),]
