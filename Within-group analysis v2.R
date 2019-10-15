##Within-group analysis script##
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


######WITHIN GROUP VARIATION########
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




#####NORMALISATION VIOLENCE#########
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

####within study####

#heatmap of intesity top varying probes
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

####Mean Correlation####14/10/19####
library(ggpubr)
betas1 <- as.data.frame(betas)

res <- cor.test(betas1$`201414140087_R01C01`, betas1$`201414140130_R01C01`, 
                method = "pearson")
res$estimate


# pdf('test.pdf')
# ggscatter(betas, x = '201414140087_R01C01', y = '201414140130_R01C01',
#           add = 'reg.line',
#           add.params = list(color = 'blue', fill = 'lightgray'),
#           conf.int = T) +
#   stat_cor(method = 'pearson') +
#   labs(title = "Correlation between UK and Indian significant sites on age (frontal)",
#        x = 'India', y = 'UK')
# dev.off()



