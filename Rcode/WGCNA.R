### Author: Pengfei Qiao (pq26@cornell.edu)
### Script to do the Weighted Gene Correlation Network Analysis in the paper

## Get scaled (sample-wise scaling, gene-wise scaling does not matter) cpm counts
library(edgeR)
counts <- read.table("counts_gene.txt")
counts <- counts[,seq(2,42,2)][,c(1:15,17:21)]
y <- DGEList(counts=counts)
y <- calcNormFactors(y)
# Remove all features that have a count of less than say 10 in more than 90% of the samples
filter_low_counts <- function (x, cutoff = 10) {
  sum(x < cutoff)
}
y.filtered <- cpm(y)[-which(apply(cpm(y), MARGIN = 1, FUN = filter_low_counts) > 18) ,]
# Variance-stabilizing transformation
library(DESeq2)
#y.vst <- varianceStabilizingTransformation(round(y.filtered))
y.log <- log(y.filtered+1, base=2) # USE this when VST doesn't provide good results
datExpr <- as.data.frame(t(y.log))

## Outlier detection
library(WGCNA)
library(cluster)
options(stringAsFactors=FALSE)

# read in biochem data
traitData <- read.csv("Cuticle_PQ_new.csv",check.names=FALSE)
traitData <- traitData[-6,]
dim(traitData)
names(traitData)
# Sum up total wax esters
traitData$totalwaxester <- apply(traitData[,20:25], 1, sum)

# Order rows of traitData to match
datTraits <- traitData[order(traitData$Sample),-1]
rownames(datTraits) <- traitData[order(traitData$Sample),1]

# Sample network based on Euclidean distance
A <- adjacency(t(datExpr), type="distance") #Cluster samples usually use Euclidean distance, clustering genes use correlation
k <- as.numeric(apply(A,2,sum)) - 1
Z.k <- scale(k)

# Designate samples as outlying based on Z.k
thresholdZ.k <- -2.5
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")

# Calculate the cluster tree using flashClust
library(flashClust)
sampleTree <- flashClust(as.dist(1-A), method = "average")
traitColors <- data.frame(numbers2colors(datTraits, signed = TRUE))
dimnames(traitColors)[[2]] <- paste(names(datTraits),"C", sep="")
datColors <- data.frame(outlierC=outlierColor, traitColors)


## Co-expression modules
# Choose a set of soft thresholding powers
powers <- c(1:30)
# choose power based on SFT criterion
disableWGCNAThreads() #To get rid of the warning in the notes of the code next line
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed") #If not run the code above, Warning message: executing %dopar% sequentially: no parallel backend registered. I think it just means it is running sequentially if no parallel backend has been registered. It will issue this warning only once. 

## Manual, stepwise module detection
# Calculate the weighted adjacency matrix, power=26
A <- adjacency(datExpr, power=26, type="signed", distFnc="bicor")
# Define a dissimilarity based on the topological overlap
dissTOM <- TOMdist(A,TOMType="signed") #Time consuming. Tried with and without TOMType="signed", identical dissTOM result. I think that's because A is positive
geneTree <- flashClust(as.dist(dissTOM), method="average")
# Define modules by cutting branches
moduleLabelsManual1 <- cutreeDynamic(dendro=geneTree, distM = dissTOM, method = "hybrid",
                                     deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
# Convert labels to colors for plotting
moduleLabelsManual <- moduleLabelsManual1
moduleColorsManual <- labels2colors(moduleLabelsManual)

# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors=moduleColorsManual)
MEs <- MEList$eigengenes

# these are the traits

totalwaxester <- as.data.frame(datTraits$totalwaxester)
names(totalwaxester) <- "total wax ester"
GS.totalwaxester <- as.numeric(cor(datExpr, totalwaxester, use="p"))
GS.totalwaxesterColor <- numbers2colors(GS.totalwaxester,signed=T)
datColors <- data.frame(moduleColorsManual,GS.totalwaxesterColor)

# Add totalwaxester
MET <- orderMEs(cbind(MEs, totalwaxester))
# Plot the relationship among the eigengenes and the trait
png("Eigengenes.png")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab=0.8, xLabelsAngle=90)
dev.off()

# Before merging highly correlated modules
png("Dendro_manual_before_merge.png")
plotDendroAndColors(geneTree,colors=datColors,groupLabels=c("Module colors","GS.total wax ester"),dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

# Automatically merge highly correlated modules
mergingThresh <- 0.1
merge <- mergeCloseModules(datExpr, moduleColorsManual,cutHeight = mergingThresh)
moduleColorsManual3 <- merge$colors #No moduleColorsAManual2 because didn't match manual with automatic, because didn't do automatic
MEsManual <- merge$newMEs

## Relating modules to traits
# Choose a module assignment
moduleColors <- moduleColorsManual3
# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEstotalalkane <- orderMEs(MEs0)
# Remove totalalkane column in datTraits
datTraits <- datTraits[,1:41]
modTraitCor <- cor(MEstotalalkane, datTraits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nSamples)
# Graphical visualization
textMatrix <- paste(signif(modTraitCor, 2),"\n(", signif(modTraitP, 1), ")", sep="")
dim(textMatrix) <- dim(modTraitCor)

# Calculate module membership (MM) values
datKME <- signedKME(datExpr, MEstotalalkane)

# Identify genes with high GS and MM
colorOfColumn <- substring(names(datKME),4)

# Get intramodular connectivity for all the genes
TOM <- TOMsimilarityFromExpr(datExpr, power=26, TOMType = "signed") #Time consuming
probes <- names(datExpr)


#########################################################
#Making figures for paper
setwd("Paper_figures")
# PCA
counts <- read.table("../counts_gene.txt", sep="\t", header=T, row.names = 1)
counts <- counts[,c(1:31,33:42)]
counts <- scale(counts)
result <- prcomp(t(counts))
# How much variantion explianed by PCs
result$sdev^2/sum(result$sdev^2)*100
#png("PCA_PC1_2.png", width = 900, height = 480)
pdf("PCA_PC1_2.pdf", width = 11.76, height=7)
# Removing outlier
plot(result$x[,1], result$x[,2], pch=16, cex = 3, xaxt='n', yaxt='n', ann=FALSE, col = rep(c("red1","red3","orange1","orange3","yellow1","yellow2","green1","green3","turquoise1","turquoise3","blue1","blue2","purple1","purple2"),3)[-32])
title(xlab = "PC1 - Developmental stage (Section 1-7)", ylab = "PC2 - Epidermal/Internal tissue", line=1, cex.lab=2)
dev.off()
# To see top or bottom ones which are epidermal
# plot(result$x[,1], result$x[,2], col = rep(c("green","red"),21))
# colors() see all the color choices

# Cuticle genes heatmap
library(edgeR)
counts <- read.table("../counts_gene.txt")
counts <- counts[,seq(2,42,2)][,c(1:15,17:21)]
y <- DGEList(counts=counts)
y <- calcNormFactors(y)
y <- data.frame(cpm(y))
# zmlip <- list.files(path="../",pattern="zmlip_*",full.names = TRUE)[-1] #-1 to get rid of the directory that matches the name
zmlip <- c("GRMZM2G101958","GRMZM2G104847","GRMZM2G088919","GRMZM2G149636","GRMZM2G177812")
zmlip_exp <- data.frame(matrix(rep(NA,7*length(zmlip)),nrow=length(zmlip)))
colnames(zmlip_exp) <- paste("Interval",c("2-4 cm", "4-6 cm","6-8 cm","8-10 cm","10-12 cm","12-14 cm", "20-22 cm"))
# rownames(zmlip_exp) <- c("ZmABC","ZmCER","ZmFAH","ZmGPAT","ZmKCS","ZmLACS","ZmLAH","ZmLTP","ZmMAGL","ZmWS","ZmWSD")
rownames(zmlip_exp) <- c("ZmLTP1","ZmLACS2","ZmGDSL","ZmKCS1","ZmABCG11")
# for (i in zmlip) {
#   a <- as.character(read.table(i)[,1])
#   x <- apply(y[which(rownames(y) %in% a),],2,sum)
#   x2 <- c(mean(x[c(1,8,15)]),mean(x[c(2,9)]),mean(x[c(3,10,16)]),mean(x[c(4,11,17)]),mean(x[c(5,12,18)]),mean(x[c(6,13,19)]),mean(x[c(7,14,20)]))
#   zmlip_exp[which(zmlip %in% i),] <- x2
# }
for (a in zmlip) {
  x <- apply(y[which(rownames(y) %in% a),],2,sum)
  x2 <- c(mean(x[c(1,8,15)]),mean(x[c(2,9)]),mean(x[c(3,10,16)]),mean(x[c(4,11,17)]),mean(x[c(5,12,18)]),mean(x[c(6,13,19)]),mean(x[c(7,14,20)]))
  zmlip_exp[which(zmlip %in% a),] <- x2
}
zmlip_exp <- t(scale(t(zmlip_exp))) #Scale is scaling columns
library(gplots)
# pdf("Cuticle_genes_heatmap2.pdf")
pdf("Cuticle_genes_heatmap3.pdf")
heatmap.2(as.matrix(zmlip_exp), Colv = FALSE, srtCol=45,dendrogram = "none", trace="none", margin = c(10,10),col=bluered, density.info="none", key.title=NA,key.xlab="Scaled expression level")
dev.off()

# Eigengenes heatmap
# colnames(MEsManual) -> c("A","B","C"....)
# selectModules <- c("lightcyan","pink","steelblue","blue","green","lightyellow","lightgreen")
MEsManual_num <- as.matrix(MEsManual)
mes <- MEsManual_num
colnames(mes) <- substring(colnames(mes),3)
# mes <- mes[,which(colnames(mes) %in% selectModules)]
mes <- mes[c(1,8,15,2,9,3,10,16,4,11,17,5,12,18,6,13,19,7,14,20),]
# Get average
mes[1,] <- (mes[1,] + mes[2,] +mes[3,])/3
mes[4,] <- (mes[4,] + mes[5,])/2
mes[6,] <- (mes[6,] + mes[7,] +mes[8,])/3
mes[9,] <- (mes[9,] + mes[10,] +mes[11,])/3
mes[12,] <- (mes[12,] + mes[13,] +mes[14,])/3
mes[15,] <- (mes[15,] + mes[16,] +mes[17,])/3
mes[18,] <- (mes[18,] + mes[19,] +mes[20,])/3
mes <- mes[seq(1,21,3),]
rownames(mes) <- paste("Interval",c("2-4 cm", "4-6 cm","6-8 cm","8-10 cm","10-12 cm","12-14 cm", "20-22 cm"))
colnames(mes) <- paste("Module", LETTERS[1:21])
mes <- t(mes)
library(gplots)
pdf("Eigengenes_heatmap2.pdf")
heatmap.2(mes, Colv = FALSE, Rowv = FALSE, dendrogram = "none", offsetRow = 0, offsetCol = 0, srtCol=45, col=bluered, trace="none", margin = c(10,10), density.info="none", key.title=NA,key.xlab="Scaled expression level")
dev.off()

# Eigengene line plot
mes <- t(mes)
mes2 <- data.frame(matrix(NA,nrow=dim(mes)[1]*dim(mes)[2],ncol=3))
mes2[,1] <- rownames(mes)
for (i in 1:dim(mes)[2]) {
  mes2[(i*dim(mes)[1]-6):(i*dim(mes)[1]),2] <- mes[,i]
  mes2[(i*dim(mes)[1]-6):(i*dim(mes)[1]),3] <- colnames(mes)[i]
}
library(ggplot2)
pdf("Eigengenes_line2.pdf")
ggplot(data=mes2,aes(x=X1,y=X2,group=X3)) + geom_line(aes(colour=X3)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Developmental stage",y="Module eigengene value") +
  scale_colour_discrete(name="Module\nmembership")
dev.off()

# Modules correlation with trait
traitData <- read.csv("Cuticle_PQ_v3_AlogmaWax.csv",check.names=FALSE)
#mw <- read.csv("molecularweight.csv",check.names = FALSE)#molecular weight, ignore warning message
elog <- rep(c(1,2.254485258,3.41124908,4.27455434,4.89023175,5.2901725,5.69280034),each=3) #Another way to adjust for cell elongation
traitData[,-1] <- traitData[,-1]*elog
traitData <- traitData[-c(6,19:21),] #Remove stage 8, remove outlier
traitData <- traitData[,-18] #Remove zero content lipid

## A shorter version
# traitData <- traitData[,c(1,8:11,13,15:17,22:33)]
# colnames(traitData)[16] <- "HFA 16:0 16-OH"
# colnames(traitData)[17] <- "HFA 16:0 9,16-diOH"
# colnames(traitData)[18] <- "HFA 18:2 9-OH"
# colnames(traitData)[19] <- "HFA 18:1 10 (9),18-diOH"
# colnames(traitData)[20] <- "HFA 18:0 9-epoxy-18-OH"
# colnames(traitData)[21] <- "DCA 18:1 9-OH"

dim(traitData)
names(traitData)
#Order rows of traitData to match
datTraits <- traitData[order(traitData$Sample),-1]
rownames(datTraits) <- traitData[order(traitData$Sample),1]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
#Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr[c(1:6,8:13,15:19),], moduleColors)$eigengenes #Remove stage 8
MEstotalalkane <- orderMEs(MEs0)
#Remove totalalkane column in datTraits
#datTraits <- datTraits[,1:41]
modTraitCor <- cor(MEstotalalkane, datTraits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nSamples)
#Graphical visualization
textMatrix <- paste(signif(modTraitCor, 2),"\n(", signif(modTraitP, 1), ")", sep="")
dim(textMatrix) <- dim(modTraitCor)

png("Cor_heatmap_merge_new_biochem.png", width = 15, height = 8, units = 'in', res = 1200) #Change this size and make readable figures!! Although the size of the file will be huge!
par(mar=c(12,14,3,3), cex=0.8)
names(MEstotalalkane) <- paste("Module",LETTERS[1:21])
#Heatmap for correlation values
labeledHeatmap(Matrix=modTraitCor, xLabels = names(datTraits), yLabels = names(MEstotalalkane),
               ySymbols = names(MEstotalalkane), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = NULL, setStdMargins = FALSE, cex.text=0.5, zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

# Hairball with c("steelblue","pink","lightcyan","lightyellow") modules
probes <- names(datExpr)
modules = c("steelblue","pink","lightcyan","lightyellow")
# Select module probes
inModule=is.finite(match(moduleColors,modules)) #Otherwise say too long...
modProbes=probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
                               nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),
                               weighted = TRUE, threshold = 0.19,nodeNames=modProbes,
                               nodeAttr = moduleColors[inModule]) #threshold was 0.02

# Hairball with blue module for phytochrome
phytochrome <- c("GRMZM2G157727","GRMZM2G181028","GRMZM2G124532","GRMZM2G092174","GRMZM2G057935","GRMZM2G129889")
probes <- names(datExpr)
modules = c("blue")
# Select module probes
inModule=is.finite(match(moduleColors,modules)) #Otherwise say too long...
modProbes=probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
                               nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),
                               weighted = TRUE, threshold = 0.2,nodeNames=modProbes,
                               nodeAttr = moduleColors[inModule]) #threshold was 0.02
