##Within-group analysis script##

load("/mnt/data1/aisha/dnam_qc/DNAm-Batch-Effect/Batch_effect_data.RData")
##libraries
library(cvequality)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(tidyr)
library(MASS)
library(tidyverse)
library(sva)
library(limma)
library(factoextra)
library(wateRmelon)
library(gplots)
library(factoextra)
library(ggpubr)
library(R.utils)

##functions
#function allows multiple plots to a page
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

#data cleaning
#there are issues with chip locations - need to redo those from the basename
QCmetrics <- separate(data = QCmetrics, col = Basename, 
                      into = c("Chip", "Chip_Position"), sep="_")
QCmetrics$Basename <- rownames(QCmetrics)

#change plate 6 of bristol to plate6b as it can overlap with manchester]
QCmetrics_br <- QCmetrics[which(QCmetrics$Institute == 'Bristol'),]
QCmetrics_no_br <- QCmetrics[which(QCmetrics$Institute != 'Bristol'),]
QCmetrics_no_br$Plate[QCmetrics_no_br$Plate == "Plate6"] <- "Plate6m" #m stands for plate 6 of manchester samples
#QCmetrics_br$Plate[QCmetrics_br$Plate == "Plate6"] <- "Plate6"
QCmetrics <- rbind(QCmetrics_no_br, QCmetrics_br)

#change the beta data to fit the the samples used now
betas <- betas[,colnames(betas) %in% rownames(QCmetrics)]

###########################################
######   WITHIN GROUP VARIATION    ########
###########################################

###Look into within group variation

#Provides pdf of each institute and pca of each variable
insts <- c('Bristol', 'KCL', 'Manchester','Oxford')
colourby = c('Plate', 'Chip', 'Chip_Position', 'Gender', 'Age', 'BR')
#creates 4 pdf with info on the insitutes sample variation and methylation distribution
for (i in 1:length(insts)){
  institute <- insts[i]
  print(institute)
  QCmetrics_ox <- QCmetrics[which(QCmetrics$Institute == institute),]
  #now plot the pcas
  pca_intensity <- as.data.frame(cbind(QCmetrics_ox$M.median, QCmetrics_ox$U.median))
  intensity_pca <- prcomp(pca_intensity, center = TRUE, scale. = TRUE)  
  print(summary(intensity_pca))
  pdfname = paste('beta', institute, 'samples differences.pdf', sep = " ")
  pdf(pdfname)
  plots <- list()
  for (i in 1:length(colourby)){
    colby = colourby[i]
    tmp = autoplot(intensity_pca, data = QCmetrics_ox, colour = colby)
    plots[[i]] <- tmp
    #plot(tm)
  }
  multiplot(plotlist = plots, cols = 2)
  
  #now plot the normality of the median Methylated & Unmethylated Intensities
  par(mfrow = c(2,2))
  hist(QCmetrics_ox$M.median, xlab = "Median M intensity", 
       main="Histogram of normal Median Methylated Intensities", cex.main=0.7)
  qqnorm(QCmetrics_ox$M.median,main="QQ plot of normal data(M,median)",pch=19)
  qqline(QCmetrics_ox$M.median)
  
  hist(QCmetrics_ox$U.median, xlab = "Median U intensity", 
       main="Histogram of Median Unmethylated Intensities", cex.main=0.7)
  qqnorm(QCmetrics_ox$U.median,main="QQ plot of normal data(U.median)",pch=19)
  qqline(QCmetrics_ox$U.median)
  
  #run within-group ANOVA
  q = list()
  q[[1]] = ggplot(QCmetrics_ox, aes(QCmetrics_ox$Plate, QCmetrics_ox$M.median)) +
    geom_violin(fill = "lightblue3") +
    geom_boxplot() +
    labs(x = "Plates", y = "M.median intensity ",
         title ="Methylated median intesity across plates") 
  
  q[[2]] = ggplot(QCmetrics_ox, aes(QCmetrics_ox$Plate, QCmetrics_ox$U.median)) +
    geom_violin(fill = "thistle1") +
    geom_boxplot()+
    labs(x = "Plates", y = "U.median intensity ",
         title ="Unmethylated median intesity across plates")
  multiplot(plotlist = q, cols = 1)
  dev.off()
}



###########################################
######   NORMALISATION VIOLENCE   #########
###########################################
###within plate
#use qual function on different plates -betas

insts <- c('Bristol', 'KCL', 'Manchester','Oxford')
dattype = list(betas, m_intensitites, u_intensities)

#this goes through three loops: For each institue, compare each plate to another plate, 
#because of the way this loop was written take a couple of pages will be replicates 

for (i in 1:length(insts)){
  institute <- insts[i]
  print(institute)
  QCmetrics_ox <- QCmetrics[which(QCmetrics$Institute == institute),] #here the old levels are retained
  QCmetrics_ox$Plate <- factor(QCmetrics_ox$Plate)
  plates = levels(QCmetrics_ox$Plate)
  print(plates)
  plots = list()
  pdfname = paste('Normlisation Violence', institute, '.pdf', sep = " ")
  pdf(pdfname, onefile = T)
  for (j in 1:length(plates)){
    plt = plates[j]
    for (k in 1:length(plates)){
      plt1 = plates[k]
      if (plt == plt1) next
      print(plt)
      print(plt1)
      v1 <- QCmetrics_ox[which(QCmetrics_ox$Plate == plt),"Basename"]
      v2 <- QCmetrics_ox[which(QCmetrics_ox$Plate == plt1),"Basename"]
      #v1 and v2 can differ in lengths - randomly takes samples from the longer plate to match lengths
      print(length(v1))
      print(length(v2))
      if(v1 > v2){
        v2 <- sample(v2, length(v1))
        print(length(v1))
        print(length(v2))
      } else  if (v1 < v2) {
        v1 <- sample(v1, length(v2))
        print(length(v1))
        print(length(v2))
      } else next 
      tmpdat_v1 <- betas[,colnames(betas) %in% v1]
      tmpdat_v2 <- betas[,colnames(betas) %in% v2]
      normv <- qual(tmpdat_v1, tmpdat_v2)
      plot(normv[,1:2], main = paste(plt, 'and', plt1, sep = " "))
    }
  }
  dev.off()
}



####within institute#####
#use qual function on different inst -betas
v1 <- filter(QCmetrics, QCmetrics$Institute == 'Bristol')[,8]
v2 <- filter(QCmetrics, QCmetrics$Institute == 'KCL')[1:140,8]
betas_v1 <- betas[,colnames(betas) %in% v1]
betas_v2 <- betas[,colnames(betas) %in% v2]
normv1 <- qual(betas_v1, betas_v2)
plot(normv[,1:2], main = 'Bristol vs KCL betas')
#use qual function on different inst - m_intensities
m_v1 <- m_intensities[,colnames(m_intensities) %in% v1]
m_v2 <- m_intensities[,colnames(m_intensities) %in% v2]
normv2 <- qual(m_v1, m_v2)
plot(normv[,1:2], main = 'Bristol vs KCL m_intensities')
#use qual function on different inst - u_intensities
u_v1 <- u_intensities[,colnames(u_intensities) %in% v1]
u_v2 <- u_intensities[,colnames(u_intensities) %in% v2]
normv3 <- qual(u_v1, u_v2)
plot(normv[,1:2], main = 'Bristol vs KCL u_intensities')

dat = list(normv1, normv2, normv3) #set what dataframe you want to test
dattype = c("betas", "m_intensitites", "u_intensities") #what titles are these dataframes from?
pdf('institute rmsd hist.pdf')
par(mfrow = c(2,2))
for (i in 1:3){
  tmpdat <- dat[i]
  tmpdat <- as.data.frame(tmpdat)
  hist(tmpdat$rmsd,
       main = paste(dattype[i], "rmsd between pl9 & pl10", sep =" "),
       xlab = "rmsd")
}
dev.off()



####heatmap of intesity top varying probes####
sigma<-apply(betas_ox2, 1, sd)# this is calculation the standard deviation of all probes
plot(hclust(dist(t(betas[order(sigma, decreasing = TRUE)[1:500],]))), 
     labels = paste(pheno$Sex, sep = "_"), 
     cex = 0.68, main = "") # this code is clustering the 5000 probes with the largest SD i.e the most variable probes in this data set

##heatmao of top varying probes
sigma<-apply(betas, 1, sd)
cell_rows <- as.data.frame(pheno$Gender, rownames(pheno))
colnames(cell_rows)<- "Gender"
cell_age <- as.data.frame(pheno$Age, rownames(pheno))
colnames(cell_age)<- "Age"
cell_inst <- as.data.frame(pheno$Institute, rownames(pheno))
colnames(cell_inst)<- "Institute"
cell_pl <- as.data.frame(pheno$Plate, rownames(pheno))
colnames(cell_pl)<- "Plate"
cell_rows <- cbind(cell_rows, cell_age, cell_inst, cell_pl)

library(pheatmap)
pdf('varying probes.pdf')
pheatmap(betas[order(sigma, decreasing = TRUE)[1:500],], scale = 'row', 
         annotation_col = cell_rows,
         cluster_rows = F,
         show_rownames = F,
         show_colnames = F,
         main="Top varying probes")
dev.off()

##heatmap of top varying probes, remove NA and Blank institution
pheno <- pheno[which(pheno$Institute != 'NA'),]
pheno <- pheno[which(pheno$Institute != 'Blank'),]
pheno <- pheno[which(pheno$Institute != ' '),]
cell_rows <- as.data.frame(pheno$Gender, rownames(pheno))
colnames(cell_rows)<- "Gender"
cell_age <- as.data.frame(pheno$Age, rownames(pheno))
colnames(cell_age)<- "Age"
cell_inst <- as.data.frame(pheno$Institute, rownames(pheno))
colnames(cell_inst)<- "Institute"
cell_pl <- as.data.frame(pheno$Plate, rownames(pheno))
colnames(cell_pl)<- "Plate"
cell_rows <- cbind(cell_rows, cell_age, cell_pl, cell_inst)

betas1 <- betas[,colnames(betas) %in% pheno$Basename]
sigma<-apply(betas1, 1, sd)
pdf('clustering probes.pdf')
plot(hclust(dist(t(betas1[order(sigma, decreasing = TRUE)[1:500],]))), 
     labels = paste(pheno$Institute, sep = "_"), cex = 0.68, main = "") 
dev.off()

library(ggplot2)
library(sets)
dend <- betas1[order(sigma, decreasing = TRUE)[1:500],]%>% # data
  scale %>% # Scale the data
  dist %>% # calculate a distance matrix, 
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.
plot(dend)
# let's add some color:
labels_colors(dend) <- 1:6
# Now each state has a color
labels_colors(dend) 



library(pheatmap)
pdf('varying probes2.pdf')
pheatmap(betas1[order(sigma, decreasing = TRUE)[1:500],], scale = 'rows',
         annotation_col = cell_rows,
         cluster_rows = F,
         show_rownames = F,
         show_colnames = F,
         main="Top varying probes")
dev.off()

pdf('varying probes3.pdf')
pheatmap(betas1[order(sigma, decreasing = TRUE)[1:500],], 
         annotation_col = cell_rows,
         cluster_rows = F,
         show_rownames = F,
         show_colnames = F,
         main="Top varying probes")
dev.off()




####Plotting boxplots of plates coloured by institue####16/10/19####
#This will measure the mean and interquatile range of the plates
vecname = list()
QCmetrics$Plate <- as.factor(QCmetrics$Plate)
for (i in 1:nlevels(QCmetrics$Plate)){
  vecname[i] <- paste('Plate', i, sep = "")
}
vecname <- unlist(vecname)
vecname2 <- insert(vecname, ats = 7, values = 'Plate6m')

QCmetrics3 <- QCmetrics
QCmetrics3$Plate <- factor(QCmetrics3$Plate, levels = vecname2)

pdf('boxplots of plate.pdf', onefile = T)
par(mar = c(6,7,6,7))
ggplot(QCmetrics3, aes(x = Plate, y = M.median, 
                      fill = Institute)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45))

ggplot(QCmetrics3, aes(x = Plate, y = U.median, 
                       fill = Institute)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45))
        
dev.off()
library(R.utils)




###########################################
##########    CORRELATION     #############
###########################################


#Mean Correlation####14/10/19####
#We try to find the correlation between the methylation of probes between plates
betas1 <- as.data.frame(betas)
cor_mat <- matrix(data = NA, nrow =  ncol(betas), ncol = ncol(betas))
cor_mat <- as.data.frame(cor_mat)
betas2 <- betas[,1:30]

res <- cor(betas2)
library("Hmisc")

res2 <- rcorr(as.matrix(betas2))
library(corrplot)
pdf('cor.pdf')
corrplot(res2[[1]], type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
dev.off()

#plots look the same as they are from the same plate. Pull out plates from different 
#institutes and compare then across each other
for (i in 1:nrow(cor_mat)){
  samp1 = betas1[,i] #take a row from the matrix
  for (k in 1:nrow(cor_mat)){
    samp2 = betas1[,k] #take another row
    res <- cor.test(samp1, samp2, method = "pearson")
    val <- res$estimate
    cor_mat[k,i] <- val
  }
}

res <- cor.test(betas1$`201414140087_R01C01`, betas1$`201414140130_R01C01`, 
                method = "pearson")
res$estimate


###########################################
##########    CHIP POSITION    ############
###########################################
####by chip####
QCmetrics <- QCmetrics %>% group_by(Chip_Position) %>% sample_n(6)
QCmetrics$Plate <- "Plate0"
qc7$Case.Control[qc7$Case.Control== '0'] <- "control"

#QCmetrics_org <- QCmetrics
QCmetrics_org -> QCmetrics
#Chip position
plate1 <- QCmetrics[which(QCmetrics$Plate == 'Plate1'),]
plate5 <- QCmetrics[which(QCmetrics$Plate == 'Plate5'),]
plate1 <- rbind(plate1, plate5)
plate1$Chip <- as.factor(plate1$Chip)
plate1$Chip_Position <- as.factor(plate1$Chip_Position)

QCmetrics <- pheno <- plate1

plates<-unique(QCmetrics$Plate)
QCmetrics$Position<-factor(QCmetrics$Position)
QCmetrics$Chip<-factor(QCmetrics$Chip, levels=rev(unique(QCmetrics$Chip))) #keeps the levels of the factor in current order rather than sorting numerically/alphabetically, also reverses this order as heatmaps plot bottom to top
require(gridExtra)

pdf('KCL and Bristol batch chip pos.pdf')
#compare all chips and positions as one plate
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
for(plate in plates){
  samples<-QCmetrics[which(QCmetrics$Plate == plate),]
  
  plateHeatmap <- ggplot(data=samples, aes(x=Position, y=Chip)) +
    scale_fill_gradientn(colours=colorRamps::matlab.like(100), limits=c(min(QCmetrics$U.median),max(QCmetrics$M.median))) +
    labs(x="", y="") +
    theme_minimal() + 
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust=1))
  
  plot1 <- plateHeatmap +
    ggtitle("Median Methylated Intensity") +
    geom_tile(aes(fill=M.median), colour = "white") +
    theme(legend.position = "none") +
    theme(plot.title = element_text(size = 10))
  
  plot2 <- plateHeatmap +
    ggtitle("Median Unmethylated Intensity") +
    geom_tile(aes(fill=U.median), colour = "white") +
    theme(legend.position = "none") +
    theme(plot.title = element_text(size = 10))
  
  legendplot<-plateHeatmap + 
    geom_tile(aes(fill=U.median), colour = "white") +
    labs(fill="Intensity") +
    scale_alpha_manual(values=c(1,1,1)) + 
    guides(alpha = guide_legend(override.aes = list(colour="black", pch=16)))
  
  legend<-g_legend(legendplot)
  grid.arrange(plot1, plot2, legend, ncol=3, widths=c(3/7, 3/7, 1/7), top=paste("", plate))
}
dev.off()

posi <- c('R02C01', 'R04C01', 'R06C01', 'R08C01')
plate1 <- plate1[(plate1$Chip_Position == 'R02C01' & plate1$Chip_Position == 'R04C01' & 
                    plate1$Chip_Position == 'R06C01' & plate1$Chip_Position == 'R08C01'),]

plate1 <- filter(plate1, plate1$Chip_Position == 'R02C01' & plate1$Chip_Position == 'R04C01' & 
                   plate1$Chip_Position == 'R06C01' & plate1$Chip_Position == 'R08C01')
plate1 <- plate1[(plate1$Chip_Position %in% posi), ] #take the same pis
plate5 <- plate5[(plate5$Chip_Position %in% posi), ]

plate5 <- plate5[order(plate5$Chip_Position),]
plate1 <- plate1[order(plate1$Chip_Position),]

toDelete <- seq(1, nrow(plate5), 2)
plate5 <- plate5[toDelete ,]
plate5 <- plate5[-9,]

betas_v1 <- betas[,colnames(betas) %in% rownames(plate1)]
betas_v2 <- betas[,colnames(betas) %in% rownames(plate5)]
normv1 <- qual(betas_v1, betas_v2)
plot(normv1[,1:2], main = 'Bristol vs KCL betas')


#rmsd hist plots
dat = list(tmpdat_v1, tmpdat_v2) #set what dataframe you want to test
dattype = c("betas", "m_intensitites", "u_intensities") #what titles are these dataframes from?
pdf('plate rmsd hist.pdf')
par(mfrow = c(2,2))
for (i in 1:2){
  tmpdat <- dat[i]
  tmpdat <- as.data.frame(tmpdat)
  hist(tmpdat$rmsd,
       main = paste(dattype[i], "rmsd between pl9 & pl10", sep =" "),
       xlab = "rmsd")
}
dev.off()

#########################
####Watermelon###########
#########################
test <- dmrse_row(betas)
test2 <- dmrse(betas)
test3 <- dmrse_col(betas)

library(quantmod)
load("/mnt/data1/BDR/QC/combined/BDR_Mset.rdat")
msetEPIC <- msetEPIC[,colnames(msetEPIC) %in% colnames(betas)]

plotmset_density<-function(mset, study=""){
  onetwo<-fData(mset)$DESIGN
  mat<-betas(mset)
  
  plot(density(mat[onetwo=="I",1], na.rm=T, bw=0.03), cex.main=0.8, main=paste(study, "Betas"), ylim=c(0, 5.2), xlab="")
  lines(density(mat[onetwo=="II",1], na.rm=T, bw=0.03), col="red")
  
  for(j in 2:ncol(mat)){
    lines(density(mat[onetwo=="I",j], na.rm=T, bw=0.03))
    lines(density(mat[onetwo=="II",j], na.rm=T, bw=0.03), col="red")
  }
  
  legend("topright", legend=c("Type I", "Type II"), lty=1, col=c("black", "red")) 
}
plotmset_density(msetEPIC, study = 'Batch')

pk <- matrix(data = NA, nrow = nrow(betas), ncol = 2)

onetwo<-fData(msetEPIC)$DESIGN
mat<-betas(msetEPIC)

for(j in 1:8){
  p <- findPeaks(mat[onetwo=="I",j])
  pk[j,1] <- sd(p)
  p2 <- findPeaks(mat[onetwo=="II",j])
  pk[j,2] <- sd(p2)
}

library(psych)

pk1 <- describe(pk[,1],type=2) 
pk1_se <- pk1$se

pk2 <- describe(pk[,2],type=2) 
pk2_se <- pk2$se

pkk <- as.data.frame(pk)
pdf('boxplot.pdf')
ggplot(pkk, aes(pkk[,1])) +
    geom_boxplot()
dev.off()

pdf('test.pdf')
plot(mat[onetwo=="I",1])
dev.off()

p = density(mat[onetwo=="I",1], na.rm=T, bw=0.03)
p1 <-p$y
p2 <- max(p1)
head(p)
#order p$y than select the two highest peaks. These two points should be the methylated
#unmethylated of type take. Place this in the matrix table
pk <- matrix(data = NA, nrow = ncol(mat), ncol = 2)
for(j in 1:ncol(mat)){
  one <- density(mat[onetwo=="I",j], na.rm=T, bw=0.03)
  oney <- one$y
  pk[j,1] <- max(oney)
  two <- density(mat[onetwo=="II",j], na.rm=T, bw=0.03)
  twoy <- two$y
  pk[j,2] <- max(twoy)
}


pk <- as.data.frame(pk)
boxplot(pk[,1], pk[,2])
sdpk1 <- sd(pk[,1])
sdpk2 <- sd(pk[,2])


#checking out fullymethylated probes
fpheno <- read.csv("/mnt/data1/BDR/QC/combined/FullyMethylatedControlSamples.csv", stringsAsFactors = F)
fopheno <- read.csv("/mnt/data1/BDR/QC/foetal_qc/FullyMethylatedControlSamples.csv", stringsAsFactors = F)

fpheno <- rbind(fpheno,fopheno)
lowintensitysamples<-which(fpheno$M.median < 2000 | fpheno$U.median < 2000)
Intensity<-rep("OK", nrow(fpheno))
Intensity[lowintensitysamples] <-"LowIntensity"

plotfactor<-as.factor(Intensity)
plot(fpheno$M.median, fpheno$U.median, pch = 16, xlab = "Median M intensity", 
     ylab = "Median U intensity", col = 'blue')
text(fpheno$M.median, fpheno$U.median, labels=fpheno$Date_Ran, cex= 0.7)

library(methylumi)
library(wateRmelon)
#Plot density lines
idatPath<-as.character(fpheno$iDAT_Location[1])
#force = T for samples that have been run on different arrays ie 450K
msetFM <- readEPIC(idatPath=idatPath, barcodes=fpheno$Basename, parallel = FALSE, force = T)
save(msetFM, file="Batch_msetFM.rdat")


plotmset_density<-function(mset, study=""){
  onetwo<-fData(mset)$DESIGN
  mat<-betas(mset)
  
  plot(density(mat[onetwo=="I",1], na.rm=T, bw=0.03), cex.main=0.8, 
       main=paste(study, "Betas"), ylim=c(0, 5.2), xlab="", xlim=c(0,1))
  lines(density(mat[onetwo=="II",1], na.rm=T, bw=0.03), col="red", xlim=c(0,1))
  
  for(j in 2:ncol(mat)){
    lines(density(mat[onetwo=="I",j], na.rm=T, bw=0.03))
    lines(density(mat[onetwo=="II",j], na.rm=T, bw=0.03), col="red")
  }
  
  legend("topright", legend=c("Type I", "Type II"), lty=1, col=c("black", "red")) 
}

plotmset_density(msetFM, study = 'Batch')

msetFM.dasen<-dasen(msetFM)
plotmset_density(msetFM.dasen, study = 'Batch')
betafm <- betas(msetFM)
betafm.dasen <- betas(msetFM.dasen)
boxplot(betafm[,1:2])

betafm <- as.data.frame(betafm)

ggplot(betafm, aes(betafm[1:2], betafm[,1:2])) +
  geom_boxplot()

#same as nasen but type I and type II backgrounds are equalized first.
#quantile normalizes methylated and unmethylated intensities separately, then calculates betas
pdf('FullyMethy2.pdf')
plotmset_density(msetFM, study = 'Batch')
plotmset_density(msetFM.dasen, study = 'Batch')
dev.off()

#Not normalised
onetwo<-fData(msetFM)$DESIGN
mat<-betas(msetFM)
pk <- matrix(data = NA, nrow = ncol(mat), ncol = 2)
for(j in 1:ncol(mat)){
  one <- density(mat[onetwo=="I",j], na.rm=T, bw=0.03)
  oney <- one$y
  pk[j,1] <- max(oney) #highest peak in type I
  two <- density(mat[onetwo=="II",j], na.rm=T, bw=0.03)
  twoy <- two$y
  pk[j,2] <- max(twoy) #highest peak in type II
}
pk <- as.data.frame(pk)
boxplot(pk[,1], pk[,2])
sdpk1 <- sd(pk[,1])
sdpk2 <- sd(pk[,2])

#Dasen Normalised
onetwo<-fData(msetFM.dasen)$DESIGN
mat<-betas(msetFM.dasen)
pk <- matrix(data = NA, nrow = ncol(mat), ncol = 2)
for(j in 1:ncol(mat)){
  one <- density(mat[onetwo=="I",j], na.rm=T, bw=0.03)
  oney <- one$y
  pk[j,1] <- max(oney) #highest peak in type I
  two <- density(mat[onetwo=="II",j], na.rm=T, bw=0.03)
  twoy <- two$y
  pk[j,2] <- max(twoy) #highest peak in type II
}
pk <- as.data.frame(pk)
boxplot(pk[,1], pk[,2])
sdpk1.n <- sd(pk[,1])
sdpk2.n <- sd(pk[,2])

png('fullmeth.png')
boxplot(betafm[,1:12])
dev.off()

png('fullmeth dasen.png')
boxplot(betafm.dasen[,1:12])
dev.off()

library(minfi)
RGset <- read.metharray.exp(base = idatPath, targets = fpheno, force = TRUE)
save(RGset, file="Batch_RGset.rdat")
MSet.raw<- preprocessRaw(RGset)
MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
controlStripPlot(RGset, controls="specificity I",sampNames = fpheno$Basename) #can chnage this



