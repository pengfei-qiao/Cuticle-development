### Author: Pengfei Qiao (pq26@cornell.edu)
### This script is to visualize networks using R


library(GGally)
library(network)
library(sna)
library(ggplot2)
require(RColorBrewer)
require(intergraph)

#Run under the directory where TOM is stored ===============================
#Color is which color you want
plotsinglenet <- function(color, threshold) {
  index <- which(moduleColorsManual3 %in% color)
  net <- TOM[index,index]
  net[net < threshold] = 0
  diag(net) <- 0
  discardindx <- which(apply(net,1,sum) == 0)
  net <- net[-discardindx,-discardindx]
  degree <- apply(net,1,sum)
  #hist(degree)
  net <- network(net, directed = FALSE)
  col <- colorRampPalette(c('blue', 'red'))(length(degree))[rank(degree)]
  ggnet2(net, color = col,node.size = 1,edge.size = 0.1) + 
    guides(color = FALSE, size = FALSE)
}
plotsinglenet("pink",0.02)
plotsinglenet("blue",0.2)



plotmultinet <- function(color, threshold) {
  index <- which(moduleColors %in% color)
  net <- TOM[index,index]
  net[net < threshold] = 0
  diag(net) <- 0
  discardindx <- which(apply(net,1,sum) == 0)
  net <- net[-discardindx,-discardindx]
  degree <- apply(net,1,sum)
  nodes <- names(datExpr)[index][-discardindx]
  net <- network(net, directed = FALSE)
  col <- c()
  for (i in nodes) {
    col <- c(col,moduleColors[which(names(datExpr) %in% i)])
  }
  ggnet2(net, color = col,node.size = 1,edge.size = 0.1) + 
    guides(color = FALSE, size = FALSE)
}
plotmultinet(c("pink","lightyellow","lightcyan","steelblue"),0.2)




#Get layout with phytochrome labeled
phytochrome <- c("GRMZM2G157727","GRMZM2G181028","GRMZM2G124532","GRMZM2G092174","GRMZM2G057935","GRMZM2G129889")
index <- which(moduleColors %in% "blue")
net <- TOM[index,index]
net[net < 0.2] = 0
diag(net) <- 0
discardindx <- which(apply(net,1,sum) == 0)
net <- net[-discardindx,-discardindx]
degree <- apply(net,1,sum)
net <- network(net, directed = FALSE)
bluenodes <- names(datExpr)[index][-discardindx]
bluenodes[which(bluenodes %in% phytochrome)]
col <- rep("seagreen3",length(bluenodes))
col[which(bluenodes %in% phytochrome)] <- "orangered"
set.seed(123)
pdf("network_blue_phytochrome.pdf")
ggnet2(net, color = col,node.size = 1,edge.size = 0.1) + 
  guides(color = FALSE, size = FALSE)
dev.off()


#Run in this directory===================================================
#Directed graph
full <- read.csv("~/Desktop/Cuticle try/CSI/FullOutput/csi_marginal.csv",row.names = 1)
known <- as.character(read.table("~/Desktop/Cuticle try/known_genes_in_GRN2.txt")[,1])
full[full < 0.025] = 0
#discardindx <- which(apply(full,1,sum) == 0)
#discardindx2 <- which(apply(full,2,sum) == 0)
#for (i in discardindx2) {
#  if (sum(discardindx == i) ==0){
#    discardindx <- c(discardindx,i)
#    }
#  }
#full <- full[-discardindx,-discardindx]

#Only subset ones with indegree to known cuticle genes
subsetindx <- which(names(full) %in% known)
full <- full[apply(full[subsetindx,],2,sum) > 0, apply(full[subsetindx,],2,sum) > 0]

degree <- apply(full,2,sum) #Calculating outdegree
net <- network(t(full),directed=TRUE) #Transpose full to make sure direction is right. It seems in original CSV column is parent, row is child, in this function is the other way, tested by below
#full[2,6] #Gene 6 not regulating gene 2
#full[6,2] #Gene 2 regulating gene 6
#col <- c("grey","red","grey","grey","grey","green",rep("grey",97)) #Should be a line from red to green. If not transpose full, then there is no arrow from red to green
col <- colorRampPalette(c('blue', 'red'))(length(degree))[rank(degree)]
col[which(names(full) %in% known)] <- "yellow"
set.seed(123)
pdf("GRN2.pdf")
ggnet2(net,node.size = 6,color = col, edge.size = 0.1,arrow.size =3,arrow.gap=0.025)
dev.off()

