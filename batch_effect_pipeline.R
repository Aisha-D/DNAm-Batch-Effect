#### Batch effect ppl ####
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

#Read in the ready made idats objects
setwd("/mnt/data1/BDR/QC/combined")
load("BDR_Mset.rdat")
load("BDR_RGset.rdat")
betas<-betas(msetEPIC)
setwd("/mnt/data1/BDR/QC/batch_analysis")

#load the correct pheno
pheno <- read.csv("/mnt/data1/BDR/QC/combined/pheno_KCL_Bris_Ox1_Ox2_Sex_Braak.csv", header = T,  stringsAsFactors = FALSE)
pheno$Oregon<-ifelse(pheno$Brain_ID=="M1481c2"|pheno$Brain_ID=="M7342c1"|pheno$Brain_ID=="M130c1"|pheno$Brain_ID=="M1205c2", TRUE,FALSE)
pheno$Empty<-pheno$Brain_ID=="EMPTY"
pheno$Control<-pheno$Brain_ID=="Blank"
rownames(pheno) <- pheno$Basename
pheno <- pheno[order(rownames(pheno)),]
betas <- betas[,order(colnames(betas))]
identical(rownames(pheno), colnames(betas))

#form a QC metrics and samples that failed at what steps
QCmetrics<-pheno
SamplesFail<-as.logical(rep("FALSE", nrow(pheno)))
Stepsummary<-as.data.frame(matrix(ncol=0, nrow=2))
rownames(Stepsummary)<-c("Failed This Step", "Total Failed")

####Check Signal Intesities####
#the median intensities for both signals are calculated and added to the QCmetrics
m_intensities<-methylated(msetEPIC)
u_intensities<-unmethylated(msetEPIC)
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)
QCmetrics<-cbind(pheno,M.median, U.median)

#plot intensities
# coloured by chip
plotfactor<-factor(pheno$SentrixID, levels=c(unique(pheno$SentrixID), "FullyMethylated", "Oregon", "Empty"))
plotfactor[pheno$Control]<-"FullyMethylated"
plotfactor[pheno$Oregon]<-"Oregon"
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
plotfactor<-factor(pheno$Plate, levels=c(unique(pheno$Plate), "FullyMethylated", "Oregon", "Empty"))
plotfactor[pheno$Control]<-"FullyMethylated"
plotfactor[pheno$Oregon]<-"Oregon"
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

#grouping is unclear - lets see if it can cluster itself
autoplot(intensity_pca, data = QCmetrics2, colour = 'Institute', frame = T) #still a lot of overlap

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


###Look into within group variation
#Oxford samples interesting have the most variation for intensitites
QCmetrics_ox <- filter(QCmetrics, QCmetrics$Institute == 'Oxford')
pca_intensity <- as.data.frame(cbind(QCmetrics_ox$M.median, QCmetrics_ox$U.median))
intensity_pca <- prcomp(pca_intensity, center = TRUE, scale. = TRUE)                     
summary(intensity_pca)
autoplot(intensity_pca, data = QCmetrics_ox, colour = 'Plate')
autoplot(intensity_pca, data = QCmetrics_ox, colour = 'Chip')
autoplot(intensity_pca, data = QCmetrics_ox, colour = 'Position')
autoplot(intensity_pca, data = QCmetrics_ox, colour = 'Gender')
autoplot(intensity_pca, data = QCmetrics_ox, colour = 'Age')
autoplot(intensity_pca, data = QCmetrics_ox, colour = 'BR')

library(factoextra)
intensity_pca
fviz_pca_contrib(intensity_pca, choice = "var", axes = 1)

screeplot(intensity_pca, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
#lda plot to seperate across plots?

#Normality test
QCmetrics_ox$Plate <- as.factor(QCmetrics_ox$Plate) 
QCmetrics_ox$intens_ratio <-as.integer(QCmetrics_ox$M.median/QCmetrics_ox$U.median)
hist(QCmetrics_ox$M.median, xlab = "Median M intensity", 
     main="Histogram of normal Median Methylated Intensities", cex.main=0.7)
qqnorm(QCmetrics_ox$M.median,main="QQ plot of normal data(M,median)",pch=19)
qqline(QCmetrics_ox$M.median)


hist(QCmetrics_ox$U.median, xlab = "Median U intensity", 
     main="Histogram of Median Unmethylated Intensities", cex.main=0.7)
qqnorm(QCmetrics_ox$U.median,main="QQ plot of normal data(U.median)",pch=19)
qqline(QCmetrics_ox$U.median)



#Run within-group ANOVA
QCmetrics_ox$Plate <- as.factor(QCmetrics_ox$Plate) 
QCmetrics_ox$intens_ratio <-as.integer(QCmetrics_ox$M.median/QCmetrics_ox$U.median)
Ox_mean <- tapply(QCmetrics_ox$intens_ratio, QCmetrics_ox$Plate, mean)
Ox_sd <- tapply(QCmetrics_ox$intens_ratio, QCmetrics_ox$Plate, sd)
results <- cbind(Ox_mean, Ox_sd)
rownames(results) <- rownames(QCmetrics_ox)
round(results, 2)
boxplot(QCmetrics_ox$intens_ratio~QCmetrics_ox$Plate, main = 'Intens Ratio by Plates', xlab = 'Plate', ylab = 'Intens Ratio')

ggplot(QCmetrics_ox, aes(QCmetrics_ox$Plate, QCmetrics_ox$M.median)) +
  geom_violin(fill = "lightblue3") +
  geom_boxplot() +
  labs(x = "Plates", y = "M.median intensity ",
        title ="Methylated median intesity across plates") 

ggplot(QCmetrics_ox, aes(QCmetrics_ox$Plate, QCmetrics_ox$U.median)) +
  geom_violin(fill = "thistle1") +
  geom_boxplot()+
  labs(x = "Plates", y = "U.median intensity ",
       title ="Unmethylated median intesity across plates")

QCmetrics_ox$intens_ratio <-QCmetrics_ox$M.median/QCmetrics_ox$U.median
anova1 <-aov(QCmetrics_ox$intens_ratio~QCmetrics_ox$Plate)
print(summary(anova1))
pht <- TukeyHSD(anova1)
print(pht)
pht1 <- as.data.frame(pht$`QCmetrics_ox$Plate`)
p_adj <- pht1$`p adj` <- as.factor(pht1$`p adj`)
ggplot(pht1, aes(rownames(pht1), pht1$diff, fill = p_adj)) +
  geom_col() +
  ylim(-0.05, 0.05) +
  labs(x = "Plates", y = "Mean Difference",
       title ="Mean differences between Plates") +
  theme(axis.text.x = element_text(angle = 45))

hist(QCmetrics_ox$intens_ratio, xlab = "Intensity ratio intensity", 
     main="Histogram of Median Unmethylated Intensities", cex.main=0.7)
qqnorm(QCmetrics_ox$U.median,main="QQ plot of normal data(U.median)",pch=19)
qqline(QCmetrics_ox$U.median)
       


##plotting means etc..
pairs(as.levelsMap(QCmetrics_ox$Plate), och = 19)
install.packages("cvequality")
library(cvequality)
library(ggplot2)
library(ggbeeswarm)
ggplot(QCmetrics_ox, 
       aes(Plate, 
          intens_ratio)) +
  geom_boxplot() +
  geom_quasirandom(alpha = 0.5) +
  theme_bw()
QCmetrics_ox2 <- QCmetrics_ox
QCmetrics_ox2$Plate <- as.factor(QCmetrics_ox2$Plate)
asymptotic_test <- with(QCmetrics_ox, asymptotic_test(Plate, intens_ratio))

#lda
library(MASS)
QCmetrics_ox.lda <- lda(Plate ~ M.median + U.median, data = QCmetrics_ox)
QCmetrics_ox.lda
QCmetrics_ox.lda.values <- predict(QCmetrics_ox.lda)
ldahist(QCmetrics_ox.lda.values$x[,1], g = Plate)
newdata <- data.frame(type = QCmetrics_ox$Plate, lda = QCmetrics_ox.lda.values$x) 
ggplot(newdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5)

#filter the betas probes that were from oxford sample
betas_ox <- as.data.frame(betas)
rownames(QCmetrics_ox) <- QCmetrics_ox$Basename
betas_ox2 <- betas_ox[, rownames(QCmetrics_ox)]
betas_ox3 <- as.data.frame(t(betas_ox2))
identical(rownames(betas_ox3), rownames(QCmetrics_ox))
betas_ox3[is.na(betas_ox3)] <- ""
library(tidyverse)
betas_ox4 <- drop_na(betas_ox3)
betas_pca <- prcomp(na.omit(betas_ox4))
fit <- princomp(na.omit(betas_ox4), cor = TRUE)

summary(intensity_pca)
autoplot(intensity_pca, data = pheno_intesity, colour = 'Plate')

##combat to correcct plate batch effect
plotlines(betas_ox2)
library(sva)
library(limma)
pheno1 <- QCmetrics_ox
batch = pheno1$Plate
modcombat = model.matrix(~1, data=pheno1)
combat_edata = ComBat(dat=pheno1, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots= F)
sum(is.na(pheno1))

modBatch = model.matrix(~as.factor(Basename) + as.factor(batch),data = environment(pheno1)) #The model matrix being used to fit the data
mod0Batch = model.matrix(~as.factor(batch),data=environment(pheno1)) #null model being compared when fitting the data
pValuesBatch = f.pvalue(betas_ox2,modBatch,mod0Batch)
plot(modBatch)

##wateRmelon
library(wateRmelon)
boxplot(log(m_intensities), las=2, cex.axis=0.8 )
pdf('beta oxford samples differences.pdf')
boxplot(log(betas_ox2), las=2, cex.axis=0.8 )
dev.off()
##bisulfite conversion batch effect???
