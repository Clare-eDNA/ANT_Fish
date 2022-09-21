setwd("~/Documents/FRASER-LAB/ARF-Position/DataAnalysis/ANT/ANT_Fish")

library(vcfR)
library(adegenet)
library(ape)
library(vegan)
library(car)
library(dartR)
library(mmod)
library(reshape2)
library(tidyverse)
library(hierfstat)

### Read in VCF file ###
vcf <- read.vcfR('finalANT.recode.vcf')
head(vcf)
head(getFIX(vcf))

### DEPTH INFO ###

#Extract the depth info
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
dp[1:4,1:6]

dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)

palette=rep_len(c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87"), 20)
# change the end number to however many samples you have

ggplot(dpf, aes(x=Sample, y=Depth)) + geom_boxplot(fill=palette) + theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_y_continuous(trans=scales::log2_trans(), expand = c(0,0), breaks=c(1, 10, 100, 1000, 5000),
                     minor_breaks=c(1:10, 2:10*10, 2:10*100, 2:5*1000)) +
  theme(panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6)) +
  theme(panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2))

### Missingness per sample ###

myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(vcf)
myMiss <- data.frame((levels(dpf$Sample)), myMiss)
colnames(myMiss) <- c('Sample', 'Missing')
# change last number to however many samples you have
ggplot(myMiss, aes(x=Sample, y=Missing)) + geom_col(fill=palette) + theme_bw() +
  labs(y='Missingness (%)') +theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=60,hjust=1))+scale_y_continuous(expand = c(0,0))


GBS<- vcfR2genlight(vcf)
GBS@ind.names
length(GBS@ind.names)

factorname <- factor(c("3","2","2","2","3","2","3","4","4","2","4","2","4","2","3","2","2","4","2","3")) # ANT
factorname <- factor(c("2","2","2","3","3","3","3","3","2","2","2","2","2")) # BVK
factorname <- factor(c("4","4","4","4","4","4","4","4","4","4","3","3","3","3","3","3","3","3","3","3")) # QMC
factorname <- factor(c("4","3","4","4","4","3","3","3","3","4","3","4","4","4","4","3","4","3","4","3","3","3","3","4","4","3","4","3","3","4","4","3","3","3","4","4","4")) #CHW
factorname <- factor(c("3","3","3","3","3","3","3","3","4","4","4","4","4","4","4","4","4","4")) #WGR

length(factorname)


GBS@pop <- factorname
GBS@pop
poplevels<-GBS@pop

### 

### Neighbor-joining tree ###

par(mar = c(4, 3, 2, 3))
par(mfrow=c(1,1))

# basic tree
tree <- njs(dist(as.matrix(GBS)))

plot(tree, "phylogram", cex=0.75, FALSE, font=2, node.pos=1, edge.width=2, label.offset=0.5)
tiplabels(pch=20, cex=1, col=c("#FF0000",  "#AAD500",  "#3200AC")[as.numeric(pop(GBS))])
axisPhylo()
title("Neighbour-joining tree of filtered ANT data")

###

GBS
data_pca <- glPca(GBS)
## if you get NAs detected in the vector means...
toRemove <- is.na(glMean(GBS, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
GBS2 <- GBS[, !toRemove]
glPca(GBS2)

## perform PCA
pca1 <- glPca(GBS2, nf=3) # select 3

## plot eigenvalues
par(mar = c(3, 3, 3, 3))
barplot(pca1$eig, main="Eigenvalues", col=heat.colors(length(pca1$eig)))

## basic plot
scatter(pca1, ratio=.2)

## plot showing groups
s.class(pca1$scores, fac = GBS2$pop, col=palette)
text(pca1$scores,labels=indNames(GBS2), cex=0.5)
add.scatter.eig(pca1$eig,3,2,1,pos="bottomleft", inset=.01, ratio=.18)

#eigen values for PCs
data_pca <- glPca(GBS2)

eig.val<-data_pca$eig
head(eig.val)
#percentages of variance for each PC
eig.perc <- 100*data_pca$eig/sum(data_pca$eig)
eig.perc
eigen<-data.frame(eig.val,eig.perc)
eigen
#writing file with both
#write.csv(eigen,file="eigen-summary.csv",row.names=TRUE,quote=FALSE)

##retrieve the 10 PC scores to be used for plotting later
pca_scores <- as.data.frame(data_pca$scores)
pca_scores$pop <- pop(GBS2)
#write.table(pca_scores, "10PCs.txt")

##PCA plotting
ggplot(pca_scores, aes(PC1, PC2)) +
  geom_point(size=3, aes(colour=pop), alpha = 0.7, shape=19) +
  scale_color_manual(values = c("#FF0000",  "#AAD500",  "#3200AC")) +
  theme_classic()

ggplot(pca_scores, aes(PC2, PC3)) +
  geom_point(size=3, aes(colour=pop), alpha = 0.7, shape=19) +
  scale_color_manual(values = c("#FF0000",  "#AAD500",  "#3200AC")) +
  theme_classic()

ggplot(pca_scores, aes(PC1, PC3)) +
  geom_point(size=3, aes(colour=pop), alpha = 0.7, shape=19) +
  scale_color_manual(values = c("#FF0000",  "#AAD500",  "#3200AC")) +
  theme_classic()

ggplot(pca_scores, aes(PC1, PC4)) +
  geom_point(size=3, aes(colour=pop), alpha = 0.7, shape=19) +
  scale_color_manual(values = c("#FF0000","#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")) +
  theme_classic()


## Notes
## For CHW, SNP filtering matters; less filtering can lead to less overlap of populations

library(adegenet)
test.dapc <- dapc(GBS)
pnw.dapc <- dapc(GBS,n.pca = 10, n.da = 2) #14 and 4 seems to work ok; 20 and 3 is also interesting

scatter(pnw.dapc, col = palette, cex = 2, legend = TRUE, clabel = F, posi.leg = "topright", scree.pca = TRUE, posi.pca = "bottomleft", posi.da = "topleft", cleg = 0.75)

par(mar = c(5, 5,5, 5))
compoplot(pnw.dapc,col = palette, posi = 'top', show.lab=TRUE,cleg = 0.5)

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(GBS)
dapc.results$indNames <- rownames(dapc.results)

# library(reshape2)
dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

head(dapc.results, n = 6)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = palette) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

tinytex::install_tinytex()
