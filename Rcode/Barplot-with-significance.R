#originally from "/Users/HomeFolder/Desktop/My Computer/Labs/Scanlon Lab/Cuticle/cuticle_chemistry/Phytochrome"
x <- read.csv("/Users/HomeFolder/Desktop/My Computer/Labs/Scanlon Lab/Cuticle/cuticle_chemistry/Phytochrome/phytochrome.csv",row.names = 1,check.names = FALSE)
x[,1:35] <- x[,1:35]/2 #Because Richard only divided by leaf area while should actually be two sides of the leaf
x <- x[,-c(9,16,21)] #Did not change fdr but should also test them for publication, just components with 0
pvals <- c()
for (i in 1:ncol(x)) {
  pvals <- c(pvals,t.test(x[1:3,i],x[4:6,i])$p.value)
}
fdrs <- p.adjust(pvals,method="fdr")
colnames(x)[which(fdrs <= 0.05)]
sig_pos <- which(fdrs <= 0.05)
nonsig_pos <- which(fdrs > 0.05)

components <- names(x)

x <- t(x)

wt <- x[,1:3]
wtmean <- apply(wt,1,mean)
wtse <- apply(wt,1,sd)/sqrt(3)

phy <- x[,4:6]
phymean <- apply(phy,1,mean)
physe <- apply(phy,1,sd)/sqrt(3)

data <- as.data.frame(matrix(NA,ncol=4,nrow=2*length(components)))
names(data) <- c("components","genotype","mean","se")
data$components <- rep(factor(components,levels=components),2) #To make the x axis labels not arranged alphabetically, but with original orders in the csv file
data$genotype <- c(rep("WT",length(components)), rep("phyB1 phyB2",length(components)))


data$mean <- c(wtmean,phymean)
data$se <- c(wtse,physe)

#Change those numbers, especially 12, 32, etc if change input csv components
waxdata <- data[c(1:32,42:73),]
cutindata <- data[c(33:41,74:82),]

library(ggplot2)
library(ggsignif)
library(gridExtra)

#Following function is to build a function which returns an element_text object, given angle (ie degrees) and positioning (ie one of x,y,top or right) information
#See: https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
#Also see: https://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
rotatedAxisElementText = function(angle,position='x'){
  angle     = angle[1]; 
  position  = position[1]
  positions = list(x=0,y=90,top=180,right=270)
  if(!position %in% names(positions))
    stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
  if(!is.numeric(angle))
    stop("'angle' must be numeric",call.=FALSE)
  rads  = (angle - positions[[ position ]])*pi/180
  hjust = 1
  vjust = 0.5
  element_text(angle=angle,vjust=vjust,hjust=hjust)
}

#Get y_position for geom_signif
wax_y <- c()
for (i in 1:32) {
  wax_y <- c(wax_y,(max(waxdata[i,3],waxdata[i+32,3])+max(waxdata[i,4],waxdata[i+32,4])+0.3))
}
cutin_y <- c()
for (i in 1:9) {
  cutin_y <- c(cutin_y,(max(cutindata[i,3],cutindata[i+9,3])+max(cutindata[i,4],cutindata[i+9,4])+0.03))
}
#Get xmin and xmax for geom_signif
cuticle_xmin <- sort(c(sig_pos-0.2,nonsig_pos),decreasing = FALSE)
wax_xmin <- cuticle_xmin[which(cuticle_xmin <= 32)]
cutin_xmin <- cuticle_xmin[which(cuticle_xmin > 32)] -32
cuticle_xmax <- sort(c(sig_pos+0.2,nonsig_pos),decreasing = FALSE)
wax_xmax <- cuticle_xmax[which(cuticle_xmax <= 32)]
cutin_xmax <- cuticle_xmax[which(cuticle_xmax > 32)] - 32

errorbarplot <- function(input) {
  ggplot(input, aes(x=components, y=mean, fill=genotype)) + 
    geom_bar(position = position_dodge(width=0.9), stat = "identity") +
    ylab(bquote('Amount ('*mu~'g /'~dm^-2*')'))+
    xlab("Wax Component") +
    # flip x and y axis to show xlab more clearly, now solved by rotating xlab: coord_flip() +
    scale_y_continuous(expand = c(0, 0),limits=c(0,21)) + #expand is to remove the padding and make axis starts at 0, limits is to show the significance layer (*) on top of the highest bar
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position = position_dodge(0.9)) +
    #To add significance level in display - see https://cran.r-project.org/web/packages/ggsignif/vignettes/intro.html
    #NAs are for not showing any significance layer about the bars
    #All the NAs for non significant will give warnings after plotting, but it's ok, don't worry about it
    geom_signif(y_position = wax_y,xmin=wax_xmin, xmax=wax_xmax,annotation=c("*","*",NA,NA,"*","*",NA,"*","*",NA,"*","*","*","*",NA,NA,"*","*",NA,"*",NA,"*",NA,"*","*",NA,NA,NA,NA,NA,NA,NA), tip_length=0) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),axis.text=element_text(size=15), axis.title=element_text(size=15), legend.key.height=unit(1,"line"), 
          legend.key.width =unit(2,"line"), legend.text=element_text(size=15),
          axis.text.x = rotatedAxisElementText(90,'x')) + #This line rotates x label
    guides(fill = guide_legend(title = "Genotype",title.theme = element_text(size = 15,angle = 0)))
}

pdf("PhyWTWax.pdf",width=12)
errorbarplot(waxdata)
dev.off()

errorbarplot <- function(input) {
  ggplot(input, aes(x=components, y=mean, fill=genotype)) + 
    geom_bar(position = position_dodge(width=0.9), stat = "identity") +
    ylab(bquote('Amount ('*mu~'g /mg)'))+
    xlab("Cutin Monomer") +
    # flip x and y axis: coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim=c(0, 1)) + #Limit ylimt but without throwing away data outside the limit. See https://stackoverflow.com/questions/25685185/limit-ggplot2-axes-without-removing-data-outside-limits-zoom
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position = position_dodge(0.9)) +
    geom_signif(y_position = cutin_y,xmin=cutin_xmin, xmax=cutin_xmax,annotation=c(rep("*",9)), tip_length=0) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),axis.text=element_text(size=15), axis.title=element_text(size=15), legend.key.height=unit(1,"line"), 
          legend.key.width =unit(2,"line"), legend.text=element_text(size=15),
          axis.text.x = rotatedAxisElementText(90,'x')) + #This line rotates x label
    guides(fill = guide_legend(title = "Genotype",title.theme = element_text(size = 15,angle = 0)))
}

pdf("PhyWTCutin.pdf",width=12)
errorbarplot(cutindata)
dev.off()

#Only draw the one that fall outside of ylim in previous figure
cutindata <- cutindata[c(4,9,13,18),] #Subset the ones that fall outside of ylim in previous figure
cutin_y <- cutin_y[c(4,9)]
cutin_xmin <- c(1:2) -.2
cutin_xmax <- c(1:2) +.2

errorbarplot <- function(input) {
  ggplot(input, aes(x=components, y=mean, fill=genotype)) + 
    geom_bar(position = position_dodge(width=0.9), stat = "identity") +
    ylab(bquote('Amount ('*mu~'g /mg)'))+
    xlab(NULL) +
    # flip x and y axis: coord_flip() +
    scale_y_continuous(expand = c(0, 0),limits = c(0,5.2)) +
    #coord_cartesian(ylim=c(0, 1)) + #Limit ylimt but without throwing away data outside the limit. See https://stackoverflow.com/questions/25685185/limit-ggplot2-axes-without-removing-data-outside-limits-zoom
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position = position_dodge(0.9)) +
    geom_signif(y_position = cutin_y,xmin=cutin_xmin, xmax=cutin_xmax,annotation=c(rep("*",2)), tip_length=0) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),axis.text=element_text(size=15), axis.title=element_text(size=15), legend.key.height=unit(1,"line"), 
          legend.key.width =unit(2,"line"), legend.text=element_text(size=15),
          axis.text.x = rotatedAxisElementText(90,'x')) + #This line rotates x label
    guides(fill = guide_legend(title = "Genotype",title.theme = element_text(size = 15,angle = 0)))
}
pdf("PhyWTCutin_2.pdf",width=4)
errorbarplot(cutindata)
dev.off()

#If you want to add sub categories to x label, or grouping, see: https://stackoverflow.com/questions/23207878/ggplot2-group-x-axis-discrete-values-into-subgroups
            