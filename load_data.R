#load the library & phenotype file
setwd("/mnt/data1/BDR/QC/batch_analysis")
library(methylumi)
library(wateRmelon)
require(gdata)
library(minfi)
library(ggplot2)
require(gridExtra)
library(plyr)
require(IlluminaHumanMethylationEPICmanifest)
library(dplyr)
library(tidyr)

#Read in the ready made idats objects
setwd("/mnt/data1/BDR/QC/combined")
load("BDR_Mset.rdat")
load("BDR_RGset.rdat")
betas<-betas(msetEPIC)
setwd("/mnt/data1/aisha/dnam_qc/DNAm-Batch-Effect")

#load the correct pheno
pheno<-read.csv("/mnt/data1/BDR/QC/combined/pheno_KCL_Bris_Ox1_Ox2_Sex_Braak.csv")
#pheno[,c(6,7)]<-NULL

## make chip name full (R often changes this to scientific notation)
pheno$Basename2<-pheno$Basename
pheno<-separate(data = pheno, col = Basename2, into = c("SentrixID", "Position"), sep="_")
pheno$Empty<-pheno$Brain_ID=="EMPTY"
pheno$Control<-pheno$Sample_ID=="Blank"

rownames(pheno) <- pheno$Basename
pheno <- pheno[order(rownames(pheno)),]
betas <- betas[,order(colnames(betas))]
identical(rownames(pheno), colnames(betas))

#form a QC metrics and samples that failed at what steps
QCmetrics<-pheno

####Check Signal Intesities####
#the median intensities for both signals are calculated and added to the QCmetrics
m_intensities<-methylated(msetEPIC)
u_intensities<-unmethylated(msetEPIC)
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)
QCmetrics<-cbind(pheno,M.median, U.median)

#Remove the blank and control samples
pheno0 <- QCmetrics[which(QCmetrics$Empty != TRUE ),]
pheno01 <- pheno0[which(pheno0$Control != TRUE ),]

#Filter for frontal brain region and braak stage less than 2
pheno1 <- pheno01[which(pheno01$BR == "Prefrontal"),]
pheno2 <- pheno1[which(pheno1$BraakTangle_numeric >=0 & pheno1$BraakTangle_numeric <=2),]
rownames(pheno2) <- pheno2$Basename
betas1 <- betas[,rownames(pheno2)]
identical(rownames(pheno2), colnames(betas1))

pheno <- pheno2
QCmetrics <- pheno
betas <- betas1

#plot intensities
# coloured by chip
plotfactor<-factor(pheno$SentrixID, levels=c(unique(pheno$SentrixID), "FullyMethylated", "Empty"))
plotfactor[pheno$Control]<-"FullyMethylated"
plotfactor[pheno$Empty]<-"Empty"

par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity", main="Histogram of Median Methylated Intensities", cex.main=0.7)
hist(U.median, xlab = "Median U intensity", main="Histogram of Median Unmethylated Intensities", cex.main=0.7)
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity", col = rainbow(nlevels(plotfactor))[factor(plotfactor)], main="Scatter plot of Signal Intensities coloured by chip")
par(xpd=TRUE)
legend("topright", levels(factor(plotfactor)), col = rainbow(nlevels(plotfactor)), pch = 10, cex=0.5)


#coloured by plate
pheno$Plate=as.character(pheno$Plate)
plotfactor<-factor(pheno$Plate, levels=c(unique(pheno$Plate), "FullyMethylated", "Empty"))
plotfactor[pheno$Control]<-"FullyMethylated"
plotfactor[pheno$Empty]<-"Empty"

par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity", main="Histogram of Median Methylated Intensities", cex.main=0.7)
hist(U.median, xlab = "Median U intensity", main="Histogram of Median Unmethylated Intensities", cex.main=0.7)
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity", col = rainbow(nlevels(plotfactor))[factor(plotfactor)], main="Scatter plot of Signal Intensities coloured by plate")
par(xpd=TRUE)
legend("topright", levels(factor(plotfactor)), col = rainbow(nlevels(plotfactor)), pch = 10, cex=0.5)
 

#plot pca of intensities
pca_intensity <- as.data.frame(cbind(QCmetrics$M.median, QCmetrics$U.median))
intensity_pca <- prcomp(pca_intensity, center = TRUE, scale. = TRUE)                     
summary(intensity_pca)
library(ggfortify)

autoplot(intensity_pca, data = QCmetrics, colour = 'Institute')


#plot PCA again but without NA, blank or missing 
QCmetrics2 <- na.omit(QCmetrics)
pca_intensity <- as.data.frame(cbind(QCmetrics2$M.median, QCmetrics2$U.median))
intensity_pca <- prcomp(pca_intensity, center = TRUE, scale. = TRUE)                     
summary(intensity_pca)
autoplot(intensity_pca, data = QCmetrics2, colour = 'Institute')
autoplot(intensity_pca, data = QCmetrics2, colour = 'Institute', frame = T)

#Let's plot the original pca but by Plate
pca_intensity <- as.data.frame(cbind(QCmetrics$M.median, QCmetrics$U.median))
intensity_pca <- prcomp(pca_intensity, center = TRUE, scale. = TRUE)                     
summary(intensity_pca)
autoplot(intensity_pca, data = QCmetrics, colour = 'Plate')

#Let's plot the original pca but by Plate
#use QCmetrics2 as there are NAs
pca_intensity <- as.data.frame(cbind(QCmetrics2$M.median, QCmetrics2$U.median))
intensity_pca <- prcomp(pca_intensity, center = TRUE, scale. = TRUE)                     
summary(intensity_pca)
autoplot(intensity_pca, data = QCmetrics2, colour = 'M.median')
autoplot(intensity_pca, data = QCmetrics2, colour = 'U.median')

#move onto the next script to filter for what we want to analyse
