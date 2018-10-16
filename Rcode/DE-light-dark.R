### Author: Pengfei Qiao (pq26@cornell.edu)
### Compare light vs dark treatment in wild types and phyB double mutant

counts <- read.delim("rawcounts.txt", header=T, row.names = 1)
samples = c()
for (i in 1:length(colnames(counts))) {
  samples[i] = strsplit(strsplit(colnames(counts)[i],".txt")[[1]],".rawcounts.")[[1]][2]
}
colnames(counts) <- samples
#split gene ids off exons
for (i in 1:length(rownames(counts))) {
  counts$uniquegeneid[i] <- strsplit(rownames(counts)[i], ".exon")[[1]][1]
}
for (i in 1:length(counts$uniquegeneid)) {
  counts$uniquegeneid2[i] <- strsplit(counts$uniquegeneid[[i]], "_E")[[1]][1]
}
#merge counts with same id together
library(plyr)
counts2 <- ddply(counts, 'uniquegeneid2', numcolwise(sum))
counts3 <- counts2[-(1:5), ] #remove all not specifically aligned counts
row.names(counts3) <- counts3$uniquegeneid2
counts3 <- counts3[, -1]
counts <- counts3
write.table(counts, file = "counts.txt", sep = "\t", quote = F, col.names = NA, row.names = TRUE)



#PCA
counts <- read.table("counts.txt", sep="\t", header=T, row.names = 1)[41:58]
counts <- scale(counts)
result <- prcomp(t(counts))
result$sdev^2/sum(result$sdev^2)*100
plot(result$x[,1], result$x[,4],type="n")
text(result$x[,1], result$x[,4],41:58)
dev.off()

d <- dist(as.matrix(counts))
hc <- hclust(d)
pdf("hclust.pdf")
plot(hc)
dev.off()

plot(result$x[,1], result$xs[,2], col = rep(c("green","red"),21))
#colors() see all the color choices

library(edgeR)
counts <- read.delim('counts.txt', header=T, row.names = 1)[,41:58]
# dir.create("DE_W22")
setwd("DE_W22")

# pariwise comparisons
DE <- function(N1, N2) {
  sectionN <- counts[,c(1,2,3,10,11,12)]
  group <- factor(c("N1", "N1", "N1", "N2", "N2", "N2"))
  group <- factor(group)
  y <- DGEList(counts = sectionN, group = group)
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) > 1) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  design <- model.matrix(~0+group)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, contrast = c(1, -1)) #Compare levels(group)[1] - levels(group)[2]
  
  #plot smear plot for DE genes
  pdfname <- paste("Smear_WT", N1, "-", N2, ".pdf", sep = "")
  de <- decideTestsDGE(lrt, p.value = 0.05)
  detags <- rownames(y)[as.logical(de)]
  pdf (file=pdfname)
  plotSmear(lrt, de.tags=detags)	# put a smearplot in that PDF with DE genes highlighted
  dev.off()
  
  #output DE genes
  txtname <- paste("DE_", N1, "-", N2, ".txt", sep="")
  lrt <- topTags(lrt, n = nrow(lrt$table))
  write.table(lrt, txtname, sep="\t", quote = FALSE, col.names = NA)
  
}

DE("WT0h","phyB0h")

#Compare with GCN, GRN candidates
gcn <- as.character(read.csv("table2.csv")[,1])
grn <- as.character(read.csv("table1.csv")[,1])
de <- read.table("DE_WTlight-dark.txt")
de <- as.character(rownames(de)[de$FDR <= 0.05 & de$logFC > 0])
length(which(gcn %in% de))
length(which(de %in% grn))
a <-read.csv("table1.csv")[which(grn %in% de),]
b <-read.csv("table2.csv")[which(gcn %in% de),]
write.table(a, "overlap-grn-de.txt",  quote=FALSE, sep = "\t", col.names = NA, row.names = TRUE)
write.table(b, "overlap-gcn-de.txt",  quote=FALSE, sep = "\t", col.names = NA, row.names = TRUE)








# Two way ANOVA
counts <- read.table("counts.txt", sep="\t", header=T, row.names = 1)[c(44:49,53:58)]
library(edgeR)
y <- DGEList(counts = counts)
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
counts <- data.frame(cpm(y))
genotype <- c(rep("WT",6),rep("phy",6))
treatment <- rep(c(rep("light6",3),rep("dark6",3)),2)
library(car)
twowayanova <- function(i) {
  x <- data.frame(exp = as.numeric(counts[i,]), genotype = as.factor(genotype), treatment = as.factor(treatment))
  x.aov = aov(exp~genotype * treatment, data = x)
  #Anova(x.aov,type="III")
  p.i <- Anova(x.aov,type="III")$`Pr(>F)`
  return(c(rownames(counts)[i],p.i[1:4]))
}

#counter = 0
#for (i in 1:dim(counts)[1]){if (! is.null(twowayanova(i))) {counter = counter +1}}

genes <- data.frame(matrix(NA,nrow=1,ncol=5))
names(genes) <- c("Gene","p-intercept","p-genotype","p-treatment","p-interaction")
rowcounter = 1
for (i in 1:dim(counts)[1]) { #Run for 3.5 min
  if (! is.null(twowayanova(i))) {genes[rowcounter,] <- twowayanova(i); rowcounter = rowcounter + 1}
}
row.names(genes) <- genes$Gene
genes <- genes[,2:5]
genes$`p-interaction` <- p.adjust(genes$`p-interaction`,method = 'fdr')
names(genes) <- c("p-intercept","p-genotype","p-treatment","fdr-interaction")
write.table(genes,"./DE_W22/genotype-treatment-interaction.txt", quote=FALSE, sep = "\t", col.names = NA, row.names = TRUE)


# empirical Bayes for interaction
counts <- read.table("../counts.txt", sep="\t", header=T, row.names = 1)[c(44:49,53:58)]
library(edgeR)
y <- DGEList(counts)
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
metadata <- data.frame(genotype=c(rep("WT",6),rep("phyB",6)),treatment=rep(c(rep("light",3),rep("dark",3)),2))
metadata$group <- interaction(metadata$genotype, metadata$treatment)
mm <- model.matrix(~ 0 + group, data = metadata)
y <- voom(y,mm,plot=TRUE)
fit <- lmFit(y,mm)
contr <- makeContrasts((groupWT.light - groupWT.dark) - (groupphyB.light - groupphyB.dark),levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
write.table(tmp2[,c("logFC","AveExpr","P.Value","adj.P.Val")], "eBayes-interaction.txt", sep="\t", quote = FALSE, col.names = NA)
