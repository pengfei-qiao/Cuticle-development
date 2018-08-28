### Author: Pengfei Qiao (pq26@cornell.edu)
### Script to make barplots in the paper

setwd("/Users/HomeFolder/Desktop/My Computer/Labs/Scanlon Lab/Meeting_files/Other meetings_talks/Student_seminar_05042017")
#To keep space in colnames, check.names = FALSE
cuticle <- read.csv("Cuticle_PQ.csv", check.names = FALSE, header = TRUE)
cuticle$Sample <- rep(2:8, each=3)
cuticle$'Total free alcohols' <- apply(cuticle[,2:6], 1, sum)
cuticle$'Total alkanes' <- apply(cuticle[,7:16], 1, sum)
cuticle$'Total aldehydes' <- apply(cuticle[,17:19], 1, sum)
cuticle$'Total wax esters' <- apply(cuticle[,20:25], 1, sum)
cuticle$'Total fatty acids' <- apply(cuticle[,27:29], 1, sum)
cuticle$'Total hydroxy fatty acids' <- apply(cuticle[,30:38],1, sum)
cuticle$'Total dicarboxylicacids' <- apply(cuticle[,39:41], 1, sum)
##Plotting
library(ggplot2)
#Save as pdf could have higher resolution!!

#Get mean and se

##For only one trait
meanse_singletrait <- function(input, cols) {
  data <- data.frame(input[,cols])
  traitnum <- length(data)
  for (i in seq(1,21,3)) {
    for (j in 1:traitnum) {
      data2 <- data #Otherwise change first element into mean, then calculate sd is not using original data!
      data[i,j] <- mean(c(data[i,j], data[i+1,j], data[i+2,j]), na.rm = TRUE)
      data[i, j+traitnum] <- sd(c(data2[i,j], data2[i+1,j], data2[i+2,j]), na.rm = TRUE)/sqrt(3)
    }
  }
  newdata <- data[seq(1,21,3),]
  newdata$section <- factor(2:8)
  return(newdata)
}

##For more than one trait - more than one bar in barplot for each section - I'm so pround of myself!
meanse_multrait <- function(input, cols) {
  #To keep space in colnames, check.names = FALSE
  data <- data.frame(input[,cols],check.names=FALSE)
  traitnum <- length(data)
  for (i in seq(1,21,3)) {
    for (j in 1:traitnum) {
      data2 <- data
      data[i,j] <- mean(c(data[i,j], data[i+1,j], data[i+2,j]), na.rm = TRUE)
      data[i, j+traitnum] <- sd(c(data2[i,j], data2[i+1,j], data2[i+2,j]), na.rm = TRUE)/sqrt(3)
    }
  }
  newdata <- data[seq(1,21,3),]
  newdata$section <- factor(2:8)
  newdata2 <- data.frame(matrix(NA, nrow=7*traitnum, ncol=1))
  newdata2$section <- rep(newdata$section, traitnum)
  newdata2$component <- rep(colnames(newdata)[1:traitnum], each=7)
  newdata2$mean <- rep(NA, 7*traitnum)
  newdata2$se <- rep(NA, 7*traitnum)
  for (k in seq(1,traitnum)) {
    rowindex <- seq(7*(k-1)+1, 7*k)
    colindex <- k
    newdata2$mean[rowindex] <- newdata[1:7, colindex]
    newdata2$se[rowindex] <- newdata[1:7, k+traitnum]
  }
  newdata2 <- newdata2[,-1]
  return(newdata2)
}

##Barplots with se as error bar - only one bar per section
totalwax <- meanse_singletrait(cuticle, 26)
colnames(totalwax) <- c("total_wax", "se1","section")
png("totalwax.png", height=480,width=720,units="px")
limits <- aes(ymax = total_wax + se1, ymin = total_wax - se1)
p <- ggplot(totalwax, aes(y=total_wax, x=section))
dodge <- position_dodge(width=0.9)
p + geom_bar(position = dodge, stat = "identity") + 
  geom_errorbar(limits, position=dodge, width = 0.25) + 
  ylab(bquote('Total wax amount ('*mu~'g /'~dm^-2*')'))+
  xlab("Section") + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text=element_text(size=30), axis.title=element_text(size=30),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
        )
dev.off()

totalcutin <- meanse_singletrait(cuticle, 42)
colnames(totalcutin) <- c("total_cutin", "se1","section")
png("total_cutin.png",height=480,width=720,units="px")
limits <- aes(ymax = total_cutin + se1, ymin = total_cutin - se1)
p <- ggplot(totalcutin, aes(y=total_cutin, x=section))
dodge <- position_dodge(width=0.9)
p + geom_bar(position = dodge, stat = "identity") + 
  geom_errorbar(limits, position=dodge, width = 0.25) + 
  ylab(bquote('Total wax amount ('*mu~'g /'~dm^-2*')'))+
  xlab("Section") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text=element_text(size=30), axis.title=element_text(size=30),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  #theme(axis.text=element_text(size=30), axis.title=element_text(size=30,face="bold"))
dev.off()

##Barplots with se as error bar - multiple bars per section - I'm so proud of myself!
errorbarplot <- function(input) {
  ggplot(input, aes(x=section, y=mean, fill=component)) + 
    geom_bar(position = position_dodge(width=0.9), stat = "identity") +
    ylab(bquote('Amount ('*mu~'g /'~dm^-2*')'))+
    xlab("Section") +
    scale_y_continuous(expand = c(0, 0)) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position = position_dodge(0.9)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),axis.text=element_text(size=30), axis.title=element_text(size=30), legend.key.height=unit(2,"line"), 
          legend.key.width =unit(2,"line"), legend.text=element_text(size=30)) +
    guides(fill = guide_legend(title = "Compound",title.theme = element_text(size = 20,angle = 0)))
}

png("all_waxes.png", width = 900, height = 480)
errorbarplot(meanse_multrait(cuticle, 43:46))
dev.off()

#png("all_waxes_excluding_alkanes.png",width = 900, height = 480)
#errorbarplot(meanse_multrait(cuticle, c(52,53,55,56)))
#dev.off()

png("all_cutin_monomers.png",width = 900, height = 480)
errorbarplot(meanse_multrait(cuticle, 47:49))
dev.off()


#ANOVA and Tukey HSD
#cuticle$Sample <- factor(cuticle$Sample)
#totalwax.lm <- lm(Total_wax ~ Sample, data=cuticle)
#anova(totalwax.lm)
#library(multcomp)
#summary(glht(totalwax.lm, linfct=mcp(Sample = "Tukey")))
