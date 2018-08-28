### Author: Pengfei Qiao (pq26@cornell.edu)
### Script to call differentially expressed genes in the paper

setwd("/Users/HomeFolder/Desktop/My Computer/Labs/Scanlon Lab/Cuticle/RNAseq/analysis/version3_final")
counts <- read.delim("counts.txt", header=T, row.names = 1)

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
write.table(counts3, "counts_gene.txt", sep = "\t", quote=F, col.names = NA)

#Doing DE
library(edgeR)
counts <- read.delim('counts_gene.txt', header=T, row.names = 1)
# dir.create("DEgenes")
setwd("/Users/HomeFolder/Desktop/My Computer/Labs/Scanlon Lab/Cuticle/RNAseq/analysis/version3_final/DEgenes")
DE <- function(Num) {
  sectionN <- counts[,c(2*(Num-1)-1, 2*(Num-1), 2*(Num-1)+13, 2*(Num-1)+14, 2*(Num-1)+27, 2*(Num-1)+28)]
  group <- factor(c("L2", "L1", "L2", "L1", "L2", "L1"))
  group <- factor(group)
  y <- DGEList(counts = sectionN, group = group)
  y <- calcNormFactors(y)
  
  #output Epidermal only genes
  EpiOnly <- sectionN[1,]
  for (i in 1:dim(sectionN)[1]) {
    if (sectionN[i,1] < y$samples[,3][1]*y$sample[,2][1]/1000000 & sectionN[i,2] > y$samples[,3][2]*y$sample[,2][2]/1000000 & sectionN[i,3] < y$samples[,3][3]*y$sample[,2][3]/1000000 & sectionN[i,4] > y$samples[,3][4]*y$sample[,2][4]/1000000 & sectionN[i,5] < y$samples[,3][5]*y$sample[,2][5]/1000000 & sectionN[i,6] > y$samples[,3][6]*y$sample[,2][6]/1000000) {
      EpiOnly <- rbind(EpiOnly,sectionN[i,])
    }
  }
  EpiOnly <- EpiOnly[-1,]
  EpiOnly <- EpiOnly[,c(2,4,6)]
  EpiOnly$mean <- (EpiOnly[,1] + EpiOnly[,2] + EpiOnly[,3])/3
  EpiOnlycolnames <- paste("section", Num, sep = "")
  colnames(EpiOnly)[4] <- EpiOnlycolnames
  EpiOnlyfilenames <- paste("EpiOnly_section", Num, ".txt", sep = "")
  write.table(EpiOnly, EpiOnlyfilenames, quote = F, sep = "\t", col.names = NA)
  
  keep <- rowSums(cpm(y) > 1) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  d <- cpm(y)
  filename <- paste("CPM_section",Num, ".txt", sep="")
  write.table(d, filename, sep="\t", quote = FALSE, col.names = NA)
  design <- model.matrix(~0+group)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, contrast = c(1, -1)) #Compare levels(group)[1] - levels(group)[2]
  
  #plot smear plot for DE genes
  pdfname <- paste("Smear_section", Num,".pdf", sep = "")
  de <- decideTestsDGE(lrt, p.value = 0.05)
  detags <- rownames(y)[as.logical(de)]
  pdf (file=pdfname)
  plotSmear(lrt, de.tags=detags)	# put a smearplot in that PDF with DE genes highlighted
  dev.off()
  
  #output DE genes
  txtname <- paste("DE_section", Num, ".txt", sep="")
  lrt <- topTags(lrt, n = nrow(lrt$table))
  write.table(lrt, txtname, sep="\t", quote = FALSE, col.names = NA)
  
}

DE(2)
DE(3)
DE(4)
DE(5)
DE(6)
DE(7)
DE(8)

#Get epidermal upregulated genes
Epi <- function(Num) {
  DEfilenames <- list.files(path=".", pattern = "DE_section", full.names = TRUE)
  allDE <- read.delim(DEfilenames[Num-1], header = TRUE, row.names = 1)
  Epiup <- allDE[allDE$FDR<= 0.05 & allDE$logFC > 0,]
  CPMfilenames <- list.files(path=".", pattern="CPM_section", full.names = TRUE)
  CPM <- read.delim(CPMfilenames[Num-1], header=TRUE, row.names = 1)
  CPM <- CPM[, c(2,4,6)]
  EpiupCPM <- CPM[which(rownames(CPM) %in% rownames(Epiup)),]
  EpiupCPMfilename <- paste("EpiCPM_section", Num, ".txt", sep="")
  EpiupCPMmeancolname <- paste("Section", Num, sep="")
  EpiupCPM$mean <- (EpiupCPM[,1] + EpiupCPM[,2] + EpiupCPM[,3])/3
  colnames(EpiupCPM)[4] <- EpiupCPMmeancolname
  write.table(EpiupCPM, EpiupCPMfilename, sep = "\t", quote = F, col.names = NA)
}

Epi(2)
Epi(3)
Epi(4)
Epi(5)
Epi(6)
Epi(7)
Epi(8)

#Move epidermal sepcific genes cpm to one file for clustering, "Epispecific" should be "Epi upregulated"
library(dplyr)
EpispecificCPMfilenames <- list.files(path=".", pattern = "EpiCPM_section", full.names = T)
EpiSpecCPM_inner <- read.delim(EpispecificCPMfilenames[1], header=T)
EpiSpecCPM_full <- read.delim(EpispecificCPMfilenames[1], header=T)
for (i in 2:length(files)) {
  a <- read.delim(files[i], header=T)
  output <- inner_join(output, a, by = "X")
  EpiSpecCPM_full <- full_join(EpiSpecCPM_full, a, by = "X")
}
rownames(EpiSpecCPM_inner) <- EpiSpecCPM_inner$X
EpiSpecCPM_inner <- EpiSpecCPM_inner[, c(5, 9, 13, 17, 21, 25, 29)]
rownames(EpiSpecCPM_full) <- EpiSpecCPM_full$X
EpiSpecCPM_full <- EpiSpecCPM_full[, c(5, 9, 13, 17, 21, 25, 29)]
write.table(EpiSpecCPM_inner, "EpiSpecCPM_inner.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
write.table(EpiSpecCPM_full, "EpiSpecCPM_full.txt", sep = "\t", row.names = T, col.names = NA, quote = F)

#Merge all epidermal only files
library(dplyr)
EpiOnlyCPMfilenames <- list.files(path=".", pattern = "EpiOnly_section", full.names = T)
EpiOnlyCPM_inner <- read.delim(EpiOnlyCPMfilenames[1], header=T)
EpiOnlyCPM_full <- read.delim(EpiOnlyCPMfilenames[1], header=T)
for (i in 2:length(EpiOnlyCPMfilenames)) {
  a <- read.delim(EpiOnlyCPMfilenames[i], header=T)
  EpiOnlyCPM_inner <- inner_join(EpiOnlyCPM_inner, a, by = "X")
  EpiOnlyCPM_full <- full_join(EpiOnlyCPM_full, a, by = "X")
}
rownames(EpiOnlyCPM_inner) <- EpiOnlyCPM_inner$X
EpiOnlyCPM_inner <- EpiOnlyCPM_inner[, c(5, 9, 13, 17, 21, 25, 29)]
rownames(EpiOnlyCPM_full) <- EpiOnlyCPM_full$X
EpiOnlyCPM_full <- EpiOnlyCPM_full[, c(5, 9, 13, 17, 21, 25, 29)]
write.table(EpiOnlyCPM_inner, "EpiOnlyCPM_inner.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
write.table(EpiOnlyCPM_full, "EpiOnlyCPM_full.txt", sep = "\t", row.names = T, col.names = NA, quote = F)



##Get all epidermal DE genes together
setwd("/Users/HomeFolder/Desktop/My Computer/Labs/Scanlon Lab/Cuticle/RNAseq/analysis/version3_final/DEgenes/")
files <- list.files(pattern = "EpiCPM_section", full.names=T)
library(dplyr)
output <- read.delim(files[1], header=T)
for (i in 2:length(files)) {
  a <- read.delim(files[i], header=T)
  output <- inner_join(output, a, by = "X")
}
rownames(output) <- output$X
output <- output[,-1] #Get rid of gene id as first column
output <- output[,-seq(4, 28, 4)] #Get ride of section mean as a column
write.table(output, "shared_epiDE_CPM.txt", sep = '\t', quote=F, col.names = NA)
