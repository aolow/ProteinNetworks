# Date: 02-06-2015
# Author: Aleksandra Olow
# Aim: identify profile of expression for the intertwined kinase-substrate MAPK/PI3K network

#rm(list=ls())
setwd("~/Tissue_Expression_BRAF")

# Human Protein Atlas accessed 02-06-2015 rna expression levels in cell lines and tissues
# contains Ensembl Gene ID

#rna <- read.csv("/Users/aolow/Documents/VantveerLab/Database_Project/Tissue_Expression_BRAF/rna.csv")


load("MAPK_PI3K_prots_ensembl.RData")
names(proteins2) <- c("Protein", "Gene")
load("subset_expression_pth.RData")


#Hive plot nodes - for later use
load("/Users/aolow/nodes2_labmeeting.RData")
nodes <- read.csv("/Users/aolow/Dropbox/MS1_Rcode/NEW_NODES3.csv", as.is=T)

## Exploration
rna_sbst[which.max(rna_sbst$Value),]
rna_sbst[which.min(rna_sbst$Value),]

# Exclude cell lines
tiss <- rna_sbst$Sample[45:76]
tiss <- "cerebral cortex"
rna_sbst2 <- rna_sbst[rna_sbst$Sample %in% tiss,]

# Normal Tissue RNA expression levels FPKM
library(lattice)
xyplot(rna_sbst$Value~rna_sbst$Sample|rna_sbst$Protein*rna_sbst$Abundance)

xyplot(rna_sbst2$Value~rna_sbst2$Protein|rna_sbst2$Sample, 
       scales=list(y=list(relation="free", rot=0),
                   x=list(relation="free", rot=90)),
       ylab="RNA \nExpression (FPKM)",
       xlab="")

# rec: save as 15x 10 in

# Log scaled
xyplot(log2(rna_sbst2$Value)~rna_sbst2$Protein|rna_sbst2$Sample, 
       scales=list(y=list(relation="free", rot=0),
                   x=list(relation="free", rot=90)),
    #   strip=FALSE,
       pch = 16,
       abline=list(h=c(2, 6), col="red"),
       ylab="RNA \nExpression [log2(FPKM)]",
       xlab="")

xyplot(log2(rna_sbst2$Value)~rna_sbst2$Protein|rna_sbst2$Sample, 
       scales=list(y=list(rot=0),
                   x=list(rot=90)),
       #   strip=FALSE,
       pch = 16,
       abline=list(h=c(2, 6), col="red"),
       ylab="RNA \nExpression [log2(FPKM)]",
       xlab="")

# Expand the basic network into all related proteins

load("unq_kin.RData")
load("unq_sub.RData")

proteins2 <- read.table("/Users/aolow/Documents/VantveerLab/Database_Project/Tissue_Expression_BRAF/proteins.txt", sep=",")
proteins2 <- unique(proteins2[,c(1,3)])
names(proteins2) <- c("Gene", "Protein")
proteins2 <- proteins2[-1,]

# Look at full network

load("subset_expression_big.RData")

rna_sbst_big <- merge(rna_sbst_big, proteins2, by="Gene")

tiss <- "cerebral cortex"
rna_sbst_big2 <- rna_sbst_big[rna_sbst_big$Sample %in% tiss,]


xyplot(log2(rna_sbst_big2$Value)~rna_sbst_big2$Protein, 
       scales=list(y=list(relation="free", rot=0),
                   x=list(relation="free", rot=90)),
       ylab="RNA \nExpression [log2(FPKM)]",
       strip=FALSE,
       abline=list(h=c(0, 6), col="red"),
       pch = 16,
       xlab="")

low_prots <- as.character(rna_sbst_big2[which(log2(rna_sbst_big2$Value)<=0),]$Protein)
high_prots <- as.character(rna_sbst_big2[which(log2(rna_sbst_big2$Value)>=6),]$Protein)

# save as 25x8in

#### CHRISTMAS-TREE PLOTS

mut_data <- read.csv("~DATABASE_VARS_PEPTIDES_COSMIC.csv")
data <- read.csv("~DATABASE_pep_var_summary.csv")
nodes <- read.csv("~NEW_NODES3.csv", as.is=T)[,1:3]

nodes <- within(nodes, SubPep <- paste(Substrate, Peptide, sep='-'))

data$SNP_NA[is.na(data$SNP_NA)==TRUE] <- 1
# data[is.na(data)==TRUE] <- 0
data <- transform(data, snp_avg = (((non.SNP+SNP_NA)+1)/(SNP+1)))

Data <- merge(nodes, data[,c(-1,-2)], by="SubPep", all=T)
Data[is.na(Data)==TRUE] <- 0
Data$VAR_pept_alter <- as.numeric(Data$VAR_pept_alter)

temp <- aggregate(VAR_pept_alter~ SubPep, data=Data, sum)
x <- aggregate(snp_avg~SubPep, data=Data, mean)
temp <- transform(temp, pep_tot=VAR_pept_alter+1)
temp <- merge(temp, x, all=T, by="SubPep")

temp <- merge(temp, Data[,c(1,3)], all.x=T, by="SubPep")


low_prots <- unique(as.character(rna_sbst_big2[which(log(rna_sbst_big2$Value)<=1),]$Protein))
high_prots <- unique(as.character(rna_sbst_big2[which(log(rna_sbst_big2$Value)>=4),]$Protein))


temp_D <- temp[temp$Substrate %in% proteins2$Protein,]
temp_D <- temp[temp$Substrate %in% low_prots,]
temp_D <- temp[temp$Substrate %in% high_prots,]
temp_D <- temp[temp$Substrate %in% proteins,]
temp_D <- unique(temp_D)


## CHRISTMAS TREE PLOT p3 for set temp_D

library(RColorBrewer)
library(ggplot2)

myPalette2 <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc2 <- scale_colour_gradientn(colours = c(myPalette2(52)) , limits=c(0, 51), 
                             name = "non-SNP/SNP \nratio")

myPalette1 <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = c("#4d4d4d", myPalette1(48)) , limits=c(0, 48), 
                             name = "non-SNP/SNP \nratio")


p3 <- ggplot(data=temp_D, aes(x=Substrate, y=VAR_pept_alter))+
  geom_point(aes(
    color = snp_avg),
    pch=17,
    alpha=0.8,
    size=4,
    #colour="black",
    width=0.8)+
  #bg = "#3182BD" ,
  #   position=position_jitter(width=0, height=0.1)) +
  theme_bw() +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x= element_text(vjust=0.9),
    axis.title.y= element_text(vjust=0.9),
    axis.text.x = element_text(angle = 90, hjust = 1))+
  #   scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70)) +
  #   scale_x_continuous(breaks=c(0, 10, 20, 30)) +
  ylab("# Peptide altering variants") +
  xlab("")+
  ggtitle("Peptides")+
  sc

p3





