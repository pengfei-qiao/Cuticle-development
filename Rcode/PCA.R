### Author: Pengfei Qiao (pq26@cornell.edu)
### Script to do the Principal Component Analysis in the paper

counts <- read.table("counts_gene.txt", sep="\t", header=T, row.names = 1)
counts <- scale(counts)
result <- prcomp(t(counts))
result$sdev/sum(result$sdev)*100
png("PCA_PC1_2.png", width = 900, height = 480)
plot(result$x[,1], result$x[,2], pch=16, cex = 3, col = rep(c("red1","red2","orange1","orange2","yellow1","yellow2","green1","green2","turquoise1","turquoise2","blue1","blue2","purple1","purple2"),3))
dev.off()

plot(result$x[,1], result$x[,2], col = rep(c("green","red"),21))
#colors() see all the color choices
