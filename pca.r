list.of.packages <- c("ggplot2","devtools","vqv/ggbiplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.dev.packages <- list.of.packages[!(list.of.packages %in% ("vqv/ggbiplot"))]
new.undev.packages <- list.of.packages[!(list.of.packages %in% c("ggplot2","devtools"))]
if(length(new.undev.packages)) install.packages(new.undev.packages)
if(length(new.dev.packages)) install_github("vqv/ggbiplot")

library(ggbiplot)
args <- commandArgs(trailingOnly = TRUE)
outDir=paste(args[1],"/png/PCA.png",sep="")
inDir=paste(args[1],"/csv/groupedBySF.csv",sep="")
dat=read.csv(inDir,header=TRUE)
pca.sf<-dat[,2]
pca.val<-dat[,c(5,6,7,10,16,17,19)]
pca.pca <- prcomp(pca.val,center=TRUE,scale.=TRUE)
g <- ggbiplot(pca.pca, obs.scale = 1, var.scale = 1, labels=pca.sf,ellipse = TRUE,circle=TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position="top")
ggsave(outDir, width=10, height=10, units="in", dpi=300)