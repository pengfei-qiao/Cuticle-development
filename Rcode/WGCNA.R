###貌似被我不断用不断改了细节。。。
##Get scaled (sample-wise scaling, gene-wise scaling does not matter) cpm counts
library(edgeR)
counts <- read.table("counts_gene.txt")
counts <- counts[,seq(2,42,2)][,c(1:15,17:21)]
y <- DGEList(counts=counts)
y <- calcNormFactors(y)
#Remove all features that have a count of less than say 10 in more than 90% of the samples
filter_low_counts <- function (x, cutoff = 10) {
  sum(x < cutoff)
}
y.filtered <- cpm(y)[-which(apply(cpm(y), MARGIN = 1, FUN = filter_low_counts) > 18) ,]
#Variance-stabilizing transformation
library(DESeq2)
#y.vst <- varianceStabilizingTransformation(round(y.filtered))
y.log <- log(y.filtered+1, base=2) # USE this when VST doesn't provide good results
datExpr <- as.data.frame(t(y.log))
# #Check VST's effect
# variances <- apply(datExpr,MARGIN=2,FUN=var)
# plot(variances)
# oldvariances <- apply(y.filtered,MARGIN=1,FUN=var)
# plot(oldvariances)

##WGCNA

##12.1 Outlier detection
library(WGCNA)
library(cluster)
options(stringAsFactors=FALSE)
#Read in the reads data set
#Compare with example dataset 
#femData <- read.csv("~/Downloads/MouseData/LiverFemale3600.csv")
#datExprFemale <- as.data.frame(t(femData[,-c(1:8)]))
#names(datExprFemale) <- femData$substanceBXH
#rownames(datExprFemale) <- names(femData)[-c(1:8)]

#read in biochem data
traitData <- read.csv("Cuticle_PQ_v2.csv",check.names=FALSE)
traitData <- traitData[-6,]
dim(traitData)
names(traitData)
#Sum up total alkanes
traitData$totalwaxester <- apply(traitData[,20:25], 1, sum)
#Compare with example dataset
#traitDat <- read.csv("~/Downloads/MouseData/ClinicalTraits.csv")
#allTraits <- traitDat[,c(2, 11:15, 17:30, 32:38)]
#Order the rows of traitData so that matches expression file

#Order rows of traitData to match
datTraits <- traitData[order(traitData$Sample),-1]
rownames(datTraits) <- traitData[order(traitData$Sample),1]

#Adjust for cell elongation
datTraits <- datTraits*c(7.385892375,12.29367057,9.506233329,5.785980406,5.777866429,5.663739532,2.831869766,7.385892375,12.29367057,9.506233329,5.785980406,5.777866429,5.663739532,2.831869766,7.385892375,9.506233329,5.785980406,5.777866429,5.663739532,2.831869766)

#Calculate cuticle amount change
x <- dim(datTraits)[2]
for (i in x:2) {
  datTraits[,i] <- datTraits[,i] - datTraits[,i-1]
}

#Sample network based on Euclidean distance
A <- adjacency(t(datExpr), type="distance") #Cluster samples usually use Euclidean distance, clustering genes use correlation
k <- as.numeric(apply(A,2,sum)) - 1
Z.k <- scale(k)

#Designate samples as outlying based on Z.k
thresholdZ.k <- -2.5
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")

#Calculate the cluster tree using flashClust
library(flashClust)
sampleTree <- flashClust(as.dist(1-A), method = "average")
traitColors <- data.frame(numbers2colors(datTraits, signed = TRUE))
dimnames(traitColors)[[2]] <- paste(names(datTraits),"C", sep="")
datColors <- data.frame(outlierC=outlierColor, traitColors)
png("trait.png")
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors, main = "Sample dendrogram and trait heatmap")
dev.off()
#No outlier!



##12.2 Co-expression modules
#12.2.1 Choose a set of soft thresholding powers
powers <- c(1:30)
#choose power based on SFT criterion
disableWGCNAThreads() #To get rid of the warning in the notes of the code next line
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed") #If not run the code above, Warning message: executing %dopar% sequentially: no parallel backend registered. I think it just means it is running sequentially if no parallel backend has been registered. It will issue this warning only once. 

#Plot the results
png("sft.png")
par(mfrow=c(1,2))
#SFT index as a function of different powers
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red")
#This line corresponds to using an R^2 cut-off of h
abline(h=0.85, col="red") #See FAQ Q6 for WGCNA. Couldn't get to 0.9. For signed should not be above 30.
#Mean connectivity as a function of different powers
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n", xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red")
#Choose power beta to be 26! Scale free topology criterion says pick a power where the saturation curve starts to saturate/come down.
dev.off()


# ###NOT DOING line 92-140!!!=============================================================================
# #12.2.2 Automatic module detection via dyamic tree cutting
# mergingThresh <- 0.25
# net <- blockwiseModules(datExpr,corType="bicor", maxBlockSize=12000,networkType="signed",power=26,minModuleSize=30,mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE, pamRespectsDendro=FALSE,saveTOMFileBase="cuticleTOM")
# moduleLabelsAutomatic <- net$colors
# # Convert labels to colors for plotting
# moduleColorsAutomatic <- labels2colors(moduleLabelsAutomatic)

# # A data frame with module eigengenes can be obtained as follows
# MEsAutomatic <- net$MEs

# #this are the traits
# totalwax <- as.data.frame(datTraits$Total_wax)
# totalcutin <- as.data.frame(datTraits$Total_cutin)
# totalalkane <- as.data.frame(datTraits$totalalkane)
# names(totalwax) <- "totalwax"
# names(totalcutin) <- "totalcutin"
# names(totalalkane) <- "totalalkane"
# # Next use this trait to define a gene significance variable
# GS.totalwax <- as.numeric(cor(datExpr,totalwax,use="p"))
# GS.totalcutin <- as.numeric(cor(datExpr, totalcutin, use="p"))
# GS.totalalkane <- as.numeric(cor(datExpr, totalalkane, use="p"))
# # This translates the numeric values into colors
# GS.totalwaxColor <- numbers2colors(GS.totalwax,signed=T)
# GS.totalcutinColor <- numbers2colors(GS.totalcutin,signed=T)
# GS.totalalkaneColor <- numbers2colors(GS.totalalkane,signed=T)
# blocknumber <- 1 #In the step generating net, if maxBlockSize = 5000, will genrate 3 blocks - can zoom into block 2 to see the dendrogram slightly better
# datColors <- data.frame(moduleColorsAutomatic,GS.totalwaxColor,GS.totalcutinColor, GS.totalalkaneColor)[net$blockGenes[[blocknumber]],]

# # Plot the dendrogram and the module colors underneath
# png("Dendro_automatic.png")
# plotDendroAndColors(net$dendrograms[[blocknumber]],colors=datColors,groupLabels=c("Module colors","GS.totalwax", "GS.totalcutin", "GS.totalalkane"),dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
# dev.off()

# #12.2.3 Blockwise module detection for large networks
# bwnet = blockwiseModules(datExpr,corType="bicor",
#                          maxBlockSize=5000,networkType="signed",power=26,minModuleSize=30,
#                          mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,
#                          pamRespectsDendro=FALSE,saveTOMFileBase="cuticleTOM-blockwise",verbose=0)
# #Relabel blockwise modules so that labels match net
# moduleLabelsBlockwise <- matchLabels(bwnet$colors, moduleLabelsAutomatic)
# moduleColorsBlockwise <- labels2colors(moduleLabelsBlockwise)
# mean(moduleLabelsBlockwise == moduleLabelsAutomatic)
# blockNumber <- 2
# png("Dendro_blockwise.png")
# plotDendroAndColors(bwnet$dendrograms[[blockNumber]],
#                     moduleColorsBlockwise[bwnet$blockGenes[[blockNumber]]],"Module colors",
#                     main=paste("Dendrogram and module colors in block",blockNumber),
#                     dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
# dev.off()

#=======================================================================================

#12.2.4 Manual, stepwise module detection
#Calculate the weighted adjacency matrix, power=26
A <- adjacency(datExpr, power=26, type="signed", distFnc="bicor")
#Define a dissimilarity based on the topological overlap
dissTOM <- TOMdist(A,TOMType="signed") #Time consuming. Tried with and without TOMType="signed", identical dissTOM result. I think that's because A is positive
geneTree <- flashClust(as.dist(dissTOM), method="average")
#Define modules by cutting branches
moduleLabelsManual1 <- cutreeDynamic(dendro=geneTree, distM = dissTOM, method = "hybrid",
                                     deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
#Relabel the manual modules so that labels match previous
#moduleLabelsManual2 <- matchLabels(moduleLabelsManual1, moduleLabelsAutomatic)
#Convert labels to colors for plotting
moduleLabelsManual <- moduleLabelsManual1
moduleColorsManual <- labels2colors(moduleLabelsManual)

#Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors=moduleColorsManual)
MEs <- MEList$eigengenes

#this are the traits
#totalwax <- as.data.frame(datTraits$'Total wax')
#totalcutin <- as.data.frame(datTraits$'Total cutin')
totalwaxester <- as.data.frame(datTraits$totalwaxester)
#names(totalwax) <- "total wax"
#names(totalcutin) <- "total cutin"
names(totalwaxester) <- "total wax ester"
# Next use this trait to define a gene significance variable
#GS.totalwax <- as.numeric(cor(datExpr,totalwax,use="p"))
#GS.totalcutin <- as.numeric(cor(datExpr, totalcutin, use="p"))
GS.totalwaxester <- as.numeric(cor(datExpr, totalwaxester, use="p"))
# This translates the numeric values into colors
#GS.totalwaxColor <- numbers2colors(GS.totalwax,signed=T)
#GS.totalcutinColor <- numbers2colors(GS.totalcutin,signed=T)
GS.totalwaxesterColor <- numbers2colors(GS.totalwaxester,signed=T)
datColors <- data.frame(moduleColorsManual,GS.totalwaxesterColor)

#Add totalalkane
MET <- orderMEs(cbind(MEs, totalwaxester))
#Plot the relationship among the eigengenes and the trait
png("Eigengenes.png")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab=0.8, xLabelsAngle=90)
dev.off()

#Before merging highly correlated modules
png("Dendro_manual_before_merge.png")
plotDendroAndColors(geneTree,colors=datColors,groupLabels=c("Module colors","GS.total wax ester"),dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

#Automatically merge highly correlated modules
mergingThresh <- 0.1
merge <- mergeCloseModules(datExpr, moduleColorsManual,cutHeight = mergingThresh)
moduleColorsManual3 <- merge$colors #No moduleColorsAManual2 because didn't match manual with automatic, because didn't do automatic
#Eigengens of the newly nerged modules
#=================================================
MEsManual <- merge$newMEs

#write out eigengenes
write.table(MEsManual,"Module Eigengenes.txt", sep="\t",quote=FALSE, col.names=NA)
#Eigengenes heatmap
selectModules <- c("black","orange","ivory","skyblue2","turquoise","green","lightsteelblue","antiquewhite4","lavenderblush3")
MEsManual_num <- as.matrix(MEsManual)
mes <- MEsManual_num
colnames(mes) <- substring(colnames(mes),3)
mes <- mes[,which(colnames(mes) %in% selectModules)]
mes <- mes[c(1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,14,21),]
#rownames(mes) <- c("Section2, Rep1", "Section2, Rep2","Section2, Rep3","Section3, Rep1","Section3, Rep2","Section3, Rep3","Section4, Rep1","Section4, Rep2","Section4, Rep3","Section5, Rep1","Section5, Rep2","Section5, Rep3","Section6, Rep1","Section6, Rep2","Section6, Rep3","Section7, Rep1", "Section7, Rep2","Section7, Rep3","Section8, Rep1","Section8, Rep2","Section8, Rep3")
#colnames(mes) <- c("Neg3","Neg4","Neg2","Neg1","Pos1","Pos2","Pos3","Pos5","Pos4") #之所以不按顺序是因为Colv要是NA 纵坐标就溢出了 但是Colv不是NA优惠打乱顺序 所以就按cluster之后的顺序标
#分别是green, lightsteelblue, antiquewhite4,lavenderblush3,black,orange,ivory,skyblue2 和turquoise module.
#Get average
for (i in seq(1,21,3)){
 mes[i,] <- (mes[i,]+mes[i+1,]+mes[i+2,])/3
}
mes <- mes[seq(1,21,3),]
rownames(mes) <- c("Section2","Section3","Section4","Section5","Section6","Section7","Section8")
library(gplots)
png("Eigengenes_heatmap.png")
#heatmap(mes,Rowv=NA)
heatmap.2(mes,Rowv=NULL,colV=NULL,trace="none",margin=c(10,10),key=TRUE,keysize=1,key.title="Expression level",density.info="none",labRow=rownames(mes),col=rainbow(256,start=0.5,end=0.65))
dev.off()
#Green module is the very dense one! - refer to the very dense part of brown module when mergethreshold 0.25

#Show the effect of module merging by plotting the original and merged colors below the tree
datColors <- data.frame(moduleColorsManual3, GS.totalalkaneColor)
png("Dendro_Manual_after_merge.png")
plotDendroAndColors(geneTree, colors = datColors, groupLabels=c("manual hybrid", "GS.totalalakne"), 
                     dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
dev.off()


#12.2.5 Relating modules to physiological traits
#Choose a module assignment
moduleColors <- moduleColorsManual3
#Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
#Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEstotalalkane <- orderMEs(MEs0)
#Remove totalalkane column in datTraits
datTraits <- datTraits[,1:41]
modTraitCor <- cor(MEstotalalkane, datTraits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nSamples)
#Graphical visualization
textMatrix <- paste(signif(modTraitCor, 2),"\n(", signif(modTraitP, 1), ")", sep="")
dim(textMatrix) <- dim(modTraitCor)
png("Cor_heatmap_adjusted_for_elongation.png", width = 15, height = 8, units = 'in', res = 1200) #Change this size and make readable figures!! Although the size of the file will be huge!
par(mar=c(12,14,3,3), cex=0.8)
#Heatmap for correlation values
labeledHeatmap(Matrix=modTraitCor, xLabels = names(datTraits), yLabels = names(MEstotalalkane),
               ySymbols = names(MEstotalalkane), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text=0.5, zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

#Calculate module membership (MM) values
datKME <- signedKME(datExpr, MEstotalalkane)

#Identify genes with high GS and MM
colorOfColumn <- substring(names(datKME),4)
png("Genes_highGS_MM.png")
selectModules <- c("black","skyblue2","turquoise","lightsteelblue","antiquewhite4","lavenderblush3")
par(mfrow=c(length(selectModules)/2, 2))
for (module in selectModules) {
  column <- match(module,colorOfColumn)
  restModule <- moduleColors == module
  verboseScatterplot(datKME[restModule, column], GS.totalalkane[restModule], xlab = paste("Module Membership", module, "module"),
                     ylab = "GS.totalalkane", main = paste("kME.", module, "vs. GS"), col=module, pch=16)
}
dev.off()

#Get phytochrome genes in lavenderblush3 module
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.

  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}
column <- match("blue",colorOfColumn)
restModule <- moduleColors == "blue"
browncolors <- addTrans(rep("steelblue4",length(rownames(datKME)[restModule])),50)
phytochrome <- c("GRMZM2G157727","GRMZM2G181028","GRMZM2G124532","GRMZM2G092174","GRMZM2G057935","GRMZM2G129889")
browncolors[match(phytochrome, rownames(datKME)[restModule])] = addTrans("black",255)
png("phytochromes_blue.png")
plot(datKME[restModule,column],GS.totalwaxester[restModule],col=browncolors,pch=16)
dev.off()

#Get lavenderblush3, orange module overlap with DEs
sharedDE <- read.delim("shared_DEs_stepwise and 23-78des.txt",header=F)

neg1overlap <- as.character(sharedDE[,3][sharedDE[,1] == "lavenderblush3"])
column <- match("lavenderblush3",colorOfColumn)
restModule <- moduleColors == "lavenderblush3"
browncolors <- addTrans(rep("steelblue4",length(rownames(datKME)[restModule])),100)
browncolors[match(neg1overlap, rownames(datKME)[restModule])] = addTrans("black",255)
png("DE_shared_lavenderblush3.png")
plot(datKME[restModule,column],GS.totalalkane[restModule],col=browncolors,pch=16)
dev.off()

pos2overlap <- as.character(sharedDE[,3][sharedDE[,1] == "orange"])
column <- match("orange",colorOfColumn)
restModule <- moduleColors == "orange"
browncolors <- addTrans(rep("orange",length(rownames(datKME)[restModule])),100)
browncolors[match(pos2overlap, rownames(datKME)[restModule])] = addTrans("darkred",255)
png("DE_shared_orange.png")
plot(datKME[restModule,column],GS.totalalkane[restModule],col=browncolors,pch=16)
dev.off()

#Plot out gene expressions of those genes!
library(ggplot2)
#Def function to get standard deviation etc
meanse_singletrait <- function(input, cols) {
  data <- data.frame(input[,cols])
  traitnum <- length(data)
  for (i in seq(1,21,3)) {
    for (j in 1:traitnum) {
      data[i,j] <- mean(c(data[i,j], data[i+1,j], data[i+2,j]), na.rm = TRUE)
      data[i, j+traitnum] <- sd(c(data[i,j], data[i+1,j], data[i+2,j]), na.rm = TRUE)/sqrt(3)
    }
  }
  newdata <- data[seq(1,21,3),]
  newdata$section <- factor(2:8)
  return(newdata)
}

#Def multiplot function for multiple plots on the same page
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

column <- match("lavenderblush3",colorOfColumn)
restModule <- moduleColors == "lavenderblush3" #Neg1 module is lavenderblush3
#Rank neg1overlap with high to low datKME, i.e. Module membership
neg1overlap_sorted <- neg1overlap[match(sort(datKME[restModule,column][match(neg1overlap, rownames(datKME)[restModule])],decreasing=T),datKME[restModule,column][match(neg1overlap, rownames(datKME)[restModule])])]
indexrange <- 1:length(neg1overlap_sorted)
png("DE_shared_lavenderblush3_expression.png")
pp <- list()
#不知道为什么 如果直接for loop pp里输出的永远是最后一个 所以就定义个函数来做
plot_function <- function(i){
xx <- meanse_singletrait(datExpr[c(1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,14,21),], which(colnames(datExpr) %in% neg1overlap_sorted[i]))
	colnames(xx) <- c(neg1overlap_sorted[i],"se1","section")
	limits <- aes(ymax = xx[,1] + xx[,2], ymin = xx[,1] - xx[,2])
	p <- ggplot(xx, aes(y=xx[,1], x=section))
	dodge <- position_dodge(width=0.9)
	return(p + geom_bar(position = dodge, stat = "identity") + 
  	geom_errorbar(limits, position=dodge, width = 0.25) + 
  	labs(y = neg1overlap_sorted[i])+
	theme(axis.text=element_text(size=10), axis.title=element_text(size=10)))
}
for (i in indexrange){pp[[i]] <- plot_function(i)}
multiplot(pp[[1]], pp[[2]], pp[[3]], pp[[4]], pp[[5]], pp[[6]], pp[[7]], pp[[8]], pp[[9]], pp[[10]], pp[[11]], pp[[12]],cols = 3)
dev.off()

column <- match("orange",colorOfColumn)
restModule <- moduleColors == "orange" #Pos2 module is orange
#Rank pos2overlap with high to low datKME, i.e. Module membership
pos2overlap_sorted <- pos2overlap[match(sort(datKME[restModule,column][match(pos2overlap, rownames(datKME)[restModule])],decreasing=T),datKME[restModule,column][match(pos2overlap, rownames(datKME)[restModule])])][1:4]
indexrange <- 1:length(pos2overlap_sorted)
png("DE_shared_orange_expression.png", height=480, width=160,units="px")
pp <- list()
#不知道为什么 如果直接for loop pp里输出的永远是最后一个 所以就定义个函数来做
plot_function <- function(i){
xx <- meanse_singletrait(datExpr[c(1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,14,21),], which(colnames(datExpr) %in% pos2overlap_sorted[i]))
	colnames(xx) <- c(pos2overlap_sorted[i],"se1","section")
	limits <- aes(ymax = xx[,1] + xx[,2], ymin = xx[,1] - xx[,2])
	p <- ggplot(xx, aes(y=xx[,1], x=section))
	dodge <- position_dodge(width=0.9)
	return(p + geom_bar(position = dodge, stat = "identity") + 
  	geom_errorbar(limits, position=dodge, width = 0.25) + 
  	labs(y =pos2overlap_sorted[i])+
	theme(axis.text=element_text(size=10), axis.title=element_text(size=10)))
}
for (i in indexrange){pp[[i]] <- plot_function(i)}
multiplot(pp[[1]], pp[[2]], pp[[3]], pp[[4]],cols = 1)
dev.off()

# Or something like this could make multiple plots on one pagelibrary(gridExtra)
# library(ggplot2)
# p <- list()
# for(i in 1:4){
#   p[[i]] <- plot
# }
# do.call(grid.arrange,p)



#Output datKME file, including module assignment
write.table(data.frame(datKME, moduleColorsManual3),"Module Eigengene based connectivity kME.txt", sep="\t",quote=FALSE, col.names=NA)

#Get intramodular connectivity for all the genes
#Time consuming. Tried with and without TOMType="signed", identical TOM result. Not sure why
TOM <- TOMsimilarityFromExpr(datExpr, power=26, TOMType = "signed") #Time consuming
probes <- names(datExpr)
for (i in names(table(moduleColorsManual3))) {
	inModule <- moduleColorsManual3 == i
	modProbes <- probes[inModule]
	modTOM <- TOM[inModule, inModule]
	kIN <- softConnectivity(datExpr[,modProbes], corFnc="bicor", type="signed",power=26)
	write.table(data.frame(gene=names(datExpr)[inModule],kIM=kIN), sprintf("%s_intramodular_connectivity.txt",i), quote=FALSE, sep="\t", row.names=FALSE)
}


#Visualizing the network
diag(dissTOM) <- NA
png("TOMplot.png")
TOMplot(dissim=dissTOM^26,dendro=geneTree,colors=moduleColorsManual3,main="Network heatmap plot, all genes")
#Warning messages:
#1: In plot.window(...) : "colors" is not a graphical parameter
#2: In plot.xy(xy, type, ...) : "colors" is not a graphical parameter
#3: In title(...) : "colors" is not a graphical parameter
#Took a little more than one hour
dev.off()

#png("TOMplot.png", width=10, height=10, units="in", res=1200)
#TOMplot(dissim=dissTOM^26,dendro=geneTree,colors=NULL,main="Network heatmap plot, all genes") #Took a little more than one hour
#dev.off()

##Check module qualities

#Check if module satisfies scale free topology
scalefreeRsq <- list()
for (i in names(table(moduleColorsManual3))) {
	a <- read.table(sprintf("%s intramodular connectivity.txt",i), sep="\t", header=T)
	b <- scaleFreePlot(sort(a$kIM))
	scalefreeRsq[i] <- b$scaleFreeRsquared
}
scalefreeRsq
#This is old number when merge threshold =0.25
#"antiquewhite4","black","brown","darkgrey","darkred","grey","grey60","ivory","lightsteelblue","magenta","mediumpurple2","orange","steelblue"
#0.92 0.84 0.81 0.86 0.95 0.32 0.89 0.72 0.63 0.73 0.56 0.64 0.97

#Check separability
separability <- data.frame(matrix(NA, nrow=13, ncol=13))
for (i in names(table(moduleColorsManual3))) {
	for (j in names(table(moduleColorsManual3))) {
		inModule1 <- moduleColorsManual3 == i
		datX1 <- datExpr[,inModule1]
		ADJ1 <- adjacency(datX1, power=26, type="signed", distFnc="bicor")
		inModule2 <- moduleColorsManual3 == j
		datX2 <- datExpr[,inModule2]
		ADJ2 <- adjacency(datX2, power=26, type="signed", distFnc="bicor")
		datX <- data.frame(datX1,datX2)
		ADJ= (bicor(datX))^26 #Warning generated: In bicor(datX) :  bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
		BetweenADJ=ADJ[1:sum(inModule1),c((sum(inModule1)+1):(sum(inModule1)+sum(inModule2)))]
		Density1=mean(as.dist(ADJ1))
		Density2=mean(as.dist(ADJ2))
		separabilityvalue=1-mean(BetweenADJ)/sqrt(Density1*Density2)
		separability[ which(names(table(moduleColorsManual3)) %in% i),  which(names(table(moduleColorsManual3)) %in% j)] <- separabilityvalue
	}
}
colnames(separability) <- names(table(moduleColorsManual3))
rownames(separability) <- names(table(moduleColorsManual3))
write.table(separability, "Separability.txt", sep="\t", quote=FALSE, col.names=NA)

#Check module density
moduleDensity <- list()
for (i in names(table(moduleColorsManual3))) {
	inModule <- moduleColorsManual3 == i
	datX <- datExpr[,inModule]
	ADJ <- adjacency(datX, power=26, type="signed", distFnc="bicor")
	Density=mean(as.dist(ADJ))
	moduleDensity[i] <- Density
}
moduleDensity #$antiquewhite4 0.04549064; black 0.02969078; brown 0.03446546; darkgrey 0.03338731; darkred 0.01616011; grey 7.918602e-05; grey60 0.02069758; ivory 0.02188132; lightsteelblue 0.0644797

#Compare modules with Differntially Expressed genes
DE <- read.table("overlap_DEs.txt",header=F,sep="\t")
restModule <- moduleColors == "lavenderblush3"
DE_lavenderblush3 <- intersect(as.character(DE[,1]),rownames(datKME)[restModule])
length(DE_lavenderblush3)

restModule <- moduleColors == "antiquewhite4"
DE_antiquewhite4 <- intersect(as.character(DE[,1]),rownames(datKME)[restModule])
length(DE_antiquewhite4)
write.csv(data.frame(DE_lavenderblush3),"DE_lavenderblush3.csv")


# #12.2.6 Output file for gene ontology analysis
# #Read in the probe annotation
# GeneAnnotation <- read.table("annotation.txt", fill = TRUE)
# #Match probes in the data set to annotation file
# probes <- names(datExpr)
# probes2annot <- match(probes, GeneAnnotation$V1)
# #Data frame with gene significances
# datGS.Traits <- data.frame(cor(datExpr, datTraits, use = "p"))
# names(datGS.Traits) <- paste("cor", names(datGS.Traits), sep="")
# datOutput <- data.frame(ProbeID=names(datExpr), GeneAnnotation[probes2annot,],moduleColors,datKME,datGS.Traits)
# #Output
# write.table(datOutput, "CuticleResults.csv", row.names = FALSE, sep = ",")


##12.3 Systems genetic analysis with NEO


##12.4 Visualizing the network
#12.4.1 Connectivity, TOM, and MDS plots
diag(dissTOM) <- NA

#Cytoscape

#Recalculate topological overlap
TOM <- TOMsimilarityFromExpr(datExpr, power=26)
for (i in c("steelblue","pink","cyan","blue","lightcyan","brown","green","darkgreen","lightyellow","lightgreen","white")){
modules = i
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
                               weighted = TRUE, threshold = 0.05,nodeNames=modProbes,
                               nodeAttr = moduleColors[inModule]) #threshold was 0.02
}

iaa <- as.character(read.table("zmlip_AuxIAA.txt")[,1])
table(moduleColors[which(names(datExpr) %in% iaa)])

bhlh <- as.character(read.table("zmlip_bHLH.txt")[,1])
table(moduleColors[which(names(datExpr) %in% bhlh)])

kcs <- as.character(read.table("zmlip_3-ketoacyl-CoA synthase.txt")[,1])
table(moduleColors[which(names(datExpr) %in% kcs)])

saur <- as.character(read.table("zmlip_SAUR.txt")[,1])
table(moduleColors[which(names(datExpr) %in% saur)])

myb <- as.character(read.table("zmlip_myb.txt")[,1])
table(moduleColors[which(names(datExpr) %in% myb)])

ap2 <- as.character(read.table("zmlip_AP2.txt")[,1])
table(moduleColors[which(names(datExpr) %in% ap2)])

jaz <- as.character(read.table("zmlip_jasmonate-zim-domain.txt")[,1])
table(moduleColors[which(names(datExpr) %in% jaz)])



#Thin out edges a little bit
brownedges <- read.table("CytoEdgebrown.txt", header=TRUE,sep="\t")
brownedges2 <- brownedges[brownedges$weight >= 0.2*max(brownedges$weight),]
write.table(brownedges2,"CytoEdgebrown0.2.txt",sep="\t", quote=FALSE, col.names=TRUE,row.names=FALSE)


#########################################################
#########################################################
##################################In case I need to rerun
library(edgeR)
counts <- read.table("counts_gene.txt")
counts <- counts[,seq(2,42,2)][,c(1:15,17:21)]
y <- DGEList(counts=counts)
y <- calcNormFactors(y)
#Remove all features that have a count of less than say 10 in more than 90% of the samples
filter_low_counts <- function (x, cutoff = 10) {
  sum(x < cutoff)
}
y.filtered <- cpm(y)[-which(apply(cpm(y), MARGIN = 1, FUN = filter_low_counts) > 18) ,]
#Variance-stabilizing transformation
library(DESeq2)
#y.vst <- varianceStabilizingTransformation(round(y.filtered))
y.log <- log(y.filtered+1, base=2) # USE this when VST doesn't provide good results
datExpr <- as.data.frame(t(y.log))

##12.1 Outlier detection
library(WGCNA)
library(cluster)
options(stringAsFactors=FALSE)
#Read in the reads data set
#Compare with example dataset 
#femData <- read.csv("~/Downloads/MouseData/LiverFemale3600.csv")
#datExprFemale <- as.data.frame(t(femData[,-c(1:8)]))
#names(datExprFemale) <- femData$substanceBXH
#rownames(datExprFemale) <- names(femData)[-c(1:8)]

#read in biochem data
traitData <- read.csv("Cuticle_PQ_new.csv",check.names=FALSE)
traitData <- traitData[-6,]
dim(traitData)
names(traitData)
#Sum up total alkanes
traitData$totalwaxester <- apply(traitData[,20:25], 1, sum)
#Compare with example dataset
#traitDat <- read.csv("~/Downloads/MouseData/ClinicalTraits.csv")
#allTraits <- traitDat[,c(2, 11:15, 17:30, 32:38)]
#Order the rows of traitData so that matches expression file

#Order rows of traitData to match
datTraits <- traitData[order(traitData$Sample),-1]
rownames(datTraits) <- traitData[order(traitData$Sample),1]

#Sample network based on Euclidean distance
A <- adjacency(t(datExpr), type="distance") #Cluster samples usually use Euclidean distance, clustering genes use correlation
k <- as.numeric(apply(A,2,sum)) - 1
Z.k <- scale(k)

#Designate samples as outlying based on Z.k
thresholdZ.k <- -2.5
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")

#Calculate the cluster tree using flashClust
library(flashClust)
sampleTree <- flashClust(as.dist(1-A), method = "average")
traitColors <- data.frame(numbers2colors(datTraits, signed = TRUE))
dimnames(traitColors)[[2]] <- paste(names(datTraits),"C", sep="")
datColors <- data.frame(outlierC=outlierColor, traitColors)




##12.2 Co-expression modules
#12.2.1 Choose a set of soft thresholding powers
powers <- c(1:30)
#choose power based on SFT criterion
disableWGCNAThreads() #To get rid of the warning in the notes of the code next line
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed") #If not run the code above, Warning message: executing %dopar% sequentially: no parallel backend registered. I think it just means it is running sequentially if no parallel backend has been registered. It will issue this warning only once. 

#12.2.4 Manual, stepwise module detection
#Calculate the weighted adjacency matrix, power=26
A <- adjacency(datExpr, power=26, type="signed", distFnc="bicor")
#Define a dissimilarity based on the topological overlap
dissTOM <- TOMdist(A,TOMType="signed") #Time consuming. Tried with and without TOMType="signed", identical dissTOM result. I think that's because A is positive
geneTree <- flashClust(as.dist(dissTOM), method="average")
#Define modules by cutting branches
moduleLabelsManual1 <- cutreeDynamic(dendro=geneTree, distM = dissTOM, method = "hybrid",
                                     deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
#Relabel the manual modules so that labels match previous
#moduleLabelsManual2 <- matchLabels(moduleLabelsManual1, moduleLabelsAutomatic)
#Convert labels to colors for plotting
moduleLabelsManual <- moduleLabelsManual1
moduleColorsManual <- labels2colors(moduleLabelsManual)

#Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors=moduleColorsManual)
MEs <- MEList$eigengenes

#this are the traits
#totalwax <- as.data.frame(datTraits$'Total wax')
#totalcutin <- as.data.frame(datTraits$'Total cutin')
totalwaxester <- as.data.frame(datTraits$totalwaxester)
#names(totalwax) <- "total wax"
#names(totalcutin) <- "total cutin"
names(totalwaxester) <- "total wax ester"
# Next use this trait to define a gene significance variable
#GS.totalwax <- as.numeric(cor(datExpr,totalwax,use="p"))
#GS.totalcutin <- as.numeric(cor(datExpr, totalcutin, use="p"))
GS.totalwaxester <- as.numeric(cor(datExpr, totalwaxester, use="p"))
# This translates the numeric values into colors
#GS.totalwaxColor <- numbers2colors(GS.totalwax,signed=T)
#GS.totalcutinColor <- numbers2colors(GS.totalcutin,signed=T)
GS.totalwaxesterColor <- numbers2colors(GS.totalwaxester,signed=T)
datColors <- data.frame(moduleColorsManual,GS.totalwaxesterColor)

#Add totalalkane
MET <- orderMEs(cbind(MEs, totalwaxester))
#Plot the relationship among the eigengenes and the trait
png("Eigengenes.png")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab=0.8, xLabelsAngle=90)
dev.off()

#Before merging highly correlated modules
png("Dendro_manual_before_merge.png")
plotDendroAndColors(geneTree,colors=datColors,groupLabels=c("Module colors","GS.total wax ester"),dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

#Automatically merge highly correlated modules
mergingThresh <- 0.1
merge <- mergeCloseModules(datExpr, moduleColorsManual,cutHeight = mergingThresh)
moduleColorsManual3 <- merge$colors #No moduleColorsAManual2 because didn't match manual with automatic, because didn't do automatic
MEsManual <- merge$newMEs

#12.2.5 Relating modules to physiological traits
#Choose a module assignment
moduleColors <- moduleColorsManual3
#Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
#Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEstotalalkane <- orderMEs(MEs0)
#Remove totalalkane column in datTraits
datTraits <- datTraits[,1:41]
modTraitCor <- cor(MEstotalalkane, datTraits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nSamples)
#Graphical visualization
textMatrix <- paste(signif(modTraitCor, 2),"\n(", signif(modTraitP, 1), ")", sep="")
dim(textMatrix) <- dim(modTraitCor)

#Calculate module membership (MM) values
datKME <- signedKME(datExpr, MEstotalalkane)

#Identify genes with high GS and MM
colorOfColumn <- substring(names(datKME),4)

#Get intramodular connectivity for all the genes
#Time consuming. Tried with and without TOMType="signed", identical TOM result. Not sure why
TOM <- TOMsimilarityFromExpr(datExpr, power=26, TOMType = "signed") #Time consuming
probes <- names(datExpr)

#########################################################
#########################################################
##################################In case I need to rerun


#########################################################
#Making figures for paper
setwd("Paper_figures")
#PCA
counts <- read.table("../counts_gene.txt", sep="\t", header=T, row.names = 1)
counts <- counts[,c(1:31,33:42)]
counts <- scale(counts)
result <- prcomp(t(counts))
#How much variantion explianed by PCs
result$sdev^2/sum(result$sdev^2)*100
#png("PCA_PC1_2.png", width = 900, height = 480)
pdf("PCA_PC1_2.pdf", width = 11.76, height=7)
#Removing outlier but not redo PCA
#plot(result$x[,1][-32], result$x[,2][-32], pch=16, cex = 3, col = rep(c("red1","red2","orange1","orange2","yellow1","yellow2","green1","green3","turquoise1","turquoise2","blue1","blue2","purple1","purple2"),3)[-32])
plot(result$x[,1], result$x[,2], pch=16, cex = 3, xaxt='n', yaxt='n', ann=FALSE, col = rep(c("red1","red3","orange1","orange3","yellow1","yellow2","green1","green3","turquoise1","turquoise3","blue1","blue2","purple1","purple2"),3)[-32])
title(xlab = "PC1 - Developmental stage (Section 1-7)", ylab = "PC2 - Epidermal/Internal tissue", line=1, cex.lab=2)
dev.off()
#To see top or bottom ones which are epidermal
#plot(result$x[,1], result$x[,2], col = rep(c("green","red"),21))
#colors() see all the color choices

#Cuticle genes heatmap
library(edgeR)
counts <- read.table("../counts_gene.txt")
counts <- counts[,seq(2,42,2)][,c(1:15,17:21)]
y <- DGEList(counts=counts)
y <- calcNormFactors(y)
y <- data.frame(cpm(y))
zmlip <- list.files(path="../",pattern="zmlip_*",full.names = TRUE)[-1] #-1 to get rid of the directory that matches the name
zmlip_exp <- data.frame(matrix(rep(NA,7*length(zmlip)),nrow=length(zmlip)))
colnames(zmlip_exp) <- paste("Section",1:7)
rownames(zmlip_exp) <- c("ZmABC","ZmCER","ZmFAH","ZmGPAT","ZmKCS","ZmLACS","ZmLAH","ZmLTP","ZmMAGL","ZmWS","ZmWSD")
for (i in zmlip) {
  a <- as.character(read.table(i)[,1])
  x <- apply(y[which(rownames(y) %in% a),],2,sum)
  x2 <- c(mean(x[c(1,8,15)]),mean(x[c(2,9)]),mean(x[c(3,10,16)]),mean(x[c(4,11,17)]),mean(x[c(5,12,18)]),mean(x[c(6,13,19)]),mean(x[c(7,14,20)]))
  zmlip_exp[which(zmlip %in% i),] <- x2
}
zmlip_exp <- t(scale(t(zmlip_exp))) #Scale is scaling columns
library(gplots)
pdf("Cuticle_genes_heatmap2.pdf")
#heatmap(zmlip_exp,Rowv=NA)
#heatmap.2(as.matrix(zmlip_exp),dendrogram = "none",trace="none",margin=c(10,10),key=TRUE,keysize=1,key.title="Expression level",density.info="none",labRow=rownames(zmlip_exp),col=rainbow(256,start=0.5, end=0.65))
heatmap.2(as.matrix(zmlip_exp), Colv = FALSE, srtCol=45,dendrogram = "none", trace="none", margin = c(10,10), density.info="none", key.title=NA,key.xlab="Scaled expression level")
dev.off()

#Eigengenes heatmap
#colnames(MEsManual) -> c("A","B","C"....)
#selectModules <- c("lightcyan","pink","steelblue","blue","green","lightyellow","lightgreen")
MEsManual_num <- as.matrix(MEsManual)
mes <- MEsManual_num
colnames(mes) <- substring(colnames(mes),3)
#mes <- mes[,which(colnames(mes) %in% selectModules)]
mes <- mes[c(1,8,15,2,9,3,10,16,4,11,17,5,12,18,6,13,19,7,14,20),]
#rownames(mes) <- c("Section2, Rep1", "Section2, Rep2","Section2, Rep3","Section3, Rep1","Section3, Rep2","Section3, Rep3","Section4, Rep1","Section4, Rep2","Section4, Rep3","Section5, Rep1","Section5, Rep2","Section5, Rep3","Section6, Rep1","Section6, Rep2","Section6, Rep3","Section7, Rep1", "Section7, Rep2","Section7, Rep3","Section8, Rep1","Section8, Rep2","Section8, Rep3")
#colnames(mes) <- c("Neg3","Neg4","Neg2","Neg1","Pos1","Pos2","Pos3","Pos5","Pos4") #之所以不按顺序是因为Colv要是NA 纵坐标就溢出了 但是Colv不是NA优惠打乱顺序 所以就按cluster之后的顺序标
#分别是green, lightsteelblue, antiquewhite4,lavenderblush3,black,orange,ivory,skyblue2 和turquoise module.
#Get average
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
#heatmap(mes,Rowv=NA)
heatmap.2(mes, Colv = FALSE, Rowv = FALSE, dendrogram = "none", offsetRow = 0, offsetCol = 0, srtCol=45, col=bluered, trace="none", margin = c(10,10), density.info="none", key.title=NA,key.xlab="Scaled expression level")
dev.off()

#Eigengene line plot
mes2 <- data.frame(matrix(NA,nrow=dim(mes)[1]*dim(mes)[2],ncol=3))
mes2[,1] <- rownames(mes)
for (i in 1:dim(mes)[2]) {
  mes2[(i*dim(mes)[1]-6):(i*dim(mes)[1]),2] <- mes[,i]
  mes2[(i*dim(mes)[1]-6):(i*dim(mes)[1]),3] <- colnames(mes)[i]
}
library(ggplot2)
pdf("Eigengenes_line.pdf")
ggplot(data=mes2,aes(x=X1,y=X2,group=X3)) + geom_line(aes(colour=X3)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Developmental stage",y="Module eigengene value") +
  scale_colour_discrete(name="Module\nmembership")
dev.off()

#Modules correlation with trait
traitData <- read.csv("Cuticle_PQ_v3_AlogmaWax.csv",check.names=FALSE)
#mw <- read.csv("molecularweight.csv",check.names = FALSE)#molecular weight, ignore warning message
elog <- rep(c(1,2.254485258,3.41124908,4.27455434,4.89023175,5.2901725,5.69280034),each=3) #Another way to adjust for cell elongation
traitData[,-1] <- traitData[,-1]*elog
traitData <- traitData[-6,]
dim(traitData)
names(traitData)
#Order rows of traitData to match
datTraits <- traitData[order(traitData$Sample),-1]
rownames(datTraits) <- traitData[order(traitData$Sample),1]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
#Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
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

#Hairball with c("steelblue","pink","lightcyan","lightyellow") modules
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

#Hairball with blue module for phytochrome
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
