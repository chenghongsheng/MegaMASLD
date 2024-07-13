setwd("./Pseudotime")


library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(SingleCellExperiment)
library(slingshot)
library(uwot)
library(mclust)
library(tradeSeq)
library(viridis)

###########color setting###########
color4gender<-c("blue",'pink')

pie(rep(1,4),col=brewer.pal(8,"Set3")[c(7,1,6,4)])
color4disease<-c(brewer.pal(8,"Set3")[c(7,1,6)],"#FF6347")

########## import files #######

metadata <- read.csv("./cNMF/metadata_NMF.final.csv", row.names=1)
metadata$condition<- factor(metadata$condition,levels=c("healthy",'NAFL','NASH','cirrhosis'))
metadata$NMFrisk<-factor(metadata$NMFrisk,levels=c("low",'high'))
metadata$Gender<-factor(metadata$Gender)

vsd_bgc <- read.csv("./Deseq2_analysis/vsd.bgc_ssGSEA.txt", row.names=1)
genesymbol<-vsd_bgc %>% select(symbol)

vsd_bgc <-vsd_bgc[,-c(831)]

pca <- prcomp(t(assays(sce)$count), scale. = FALSE)
rd1 <- as.data.frame(pca$x[,1:2])
#all(rownames(rd1)==rownames(metadata))
#rd1$condition<-metadata$condition

metadata$PC1<-as.data.frame(rd1)$PC1
metadata$PC2<-as.data.frame(rd1)$PC2



ggplot(metadata,aes(x=PC1,y=PC2,color=condition))+
  geom_point()+
  scale_color_manual(values=color4disease)

summary(metadata$PC1)
summary(metadata$PC2)

#omit outliers based on IQR

IQR.PC1<-summary(metadata$PC1)[5]-summary(metadata$PC1)[2]
Upper.PC1<-summary(metadata$PC1)[5]+1.5*IQR.PC1
Lower.PC1<-summary(metadata$PC1)[2]-1.5*IQR.PC1

IQR.PC2<-summary(metadata$PC2)[5]-summary(metadata$PC2)[2]
Upper.PC2<-summary(metadata$PC2)[5]+1.5*IQR.PC2
Lower.PC2<-summary(metadata$PC2)[2]-1.5*IQR.PC2

metadata_sub<-subset(metadata,metadata$PC1>Lower.PC1 & 
                       metadata$PC1<Upper.PC1 &
                       metadata$PC2>Lower.PC2 & 
                       metadata$PC2<Upper.PC2)

ggplot(metadata_sub,aes(x=PC1,y=PC2,color=condition))+
  geom_point()+
  scale_color_manual(values=color4disease)


vsd_bgc.sub<-select(vsd_bgc,rownames(metadata_sub))

genefilter <- vsd_bgc
genefilter$var<-rowVars(as.matrix(genefilter))
genefilter<-genefilter %>%
  arrange(desc(var)) %>%
  top_n(1000)

sce <- SingleCellExperiment(assays = List(counts = as.matrix(vsd_bgc.sub)))
sce <-sce[rownames(genefilter),]

pca <- prcomp(t(assays(sce)$count), scale. = FALSE)
rd1 <- as.data.frame(pca$x[,1:2])
#all(rownames(rd1)==rownames(metadata_sub))
#rd1$condition<-metadata$condition

ggplot(rd1,aes(PC1,PC2,color=metadata_sub$condition))+
  geom_point()+
  scale_color_manual(values=color4disease)


#construct UMAP
rd2 <- uwot::umap(t(assays(sce)$count))
colnames(rd2) <- c('UMAP1', 'UMAP2')


ggplot(rd2,aes(UMAP1,UMAP2,color=metadata_sub$condition))+
  geom_point()+
  scale_color_manual(values=color4disease)


reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)


#clustering by mclust
cl1 <- Mclust(rd1)$classification
ggplot(rd1,aes(PC1,PC2,color=cl1))+
  geom_point()
cl1<-factor(cl1)
colData(sce)$GMM <- cl1

cl1_UMAP <- Mclust(rd2)$classification
cl1_UMAP<-factor(cl1_UMAP,levels=c(3,2,9,1,4,8,5,6,7))
ggplot(rd2,aes(UMAP1,UMAP2,color=factor(cl1_UMAP)))+
  geom_point()
colData(sce)$GMM_UMAP <- cl1_UMAP



#clustering by kmeans
cl2 <- kmeans(rd1, centers = 2)$cluster
ggplot(rd1,aes(PC1,PC2,color=cl2))+
  geom_point()
cl2<-factor(cl2)
colData(sce)$kmeans <- cl2

cl2_UMAP <- kmeans(rd2, centers = 2)$cluster
ggplot(rd1,aes(PC1,PC2,color=cl2_UMAP))+
  geom_point()
cl2_UMAP<-factor(cl2_UMAP)
colData(sce)$kmeans_UMAP <- cl2_UMAP


#plot pseudotime
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'UMAP')
summary(sce$slingPseudotime_1)

plotcol <- col[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col = 'black')



col<-rev(viridis(100,option="plasma"))

metadata_sub$pst_GMM.PCA<-sce$slingPseudotime_1
metadata_sub$PC1<-as.data.frame(reducedDims(sce)$PCA)$PC1
metadata_sub$PC2<-as.data.frame(reducedDims(sce)$PCA)$PC2
metadata_sub$UMAP1<-as.data.frame(reducedDims(sce)$UMAP)$UMAP1
metadata_sub$UMAP2<-as.data.frame(reducedDims(sce)$UMAP)$UMAP2


#p<-ggplot(metadata_sub,aes(x=UMAP1,y=UMAP2))+
#  geom_point(aes(fill = pst_GMM.PCA), col = "grey70", shape = 21,size=2,alpha=0.6)+
#  scale_fill_gradientn(colors=col)+
#  theme_minimal()

#curves <- slingCurves(sce, as.df = TRUE)
#p + geom_path(data = curves %>% arrange(Order),
#             aes(group = Lineage)) 


pst_sub<-ggplot(metadata_sub,aes(x=UMAP1,y=UMAP2))+
  geom_point(aes(fill = pst_GMM.PCA), col = "grey70", shape = 21,size=2,alpha=0.6)+
  scale_fill_gradientn(colors=col)+
  theme_minimal()+
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage),linewidth=1,color='red') 

cond_pst_sub<-ggplot(metadata_sub,aes(x=UMAP1,y=UMAP2))+
  geom_point(aes(fill = condition), col = "grey70", shape = 21,size=2,alpha=0.8)+
  scale_fill_manual(values=color4disease)+
  theme_minimal()+
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage),linewidth=1,color='red')


dev.new()
pdf('pst.pdf')
pst_sub
dev.off()

dev.new()
pdf('pst_cond.pdf')
cond_pst_sub
dev.off()



## calculate take-off point for all genes ####

take_off.pt <- c()  # Initialize an empty vector to store change points

for(col in 1:41005){
  smooth <- predict(loess(vsd_bgc[,col] ~ vsd_bgc$pst, data = vsd_bgc), seq(0,7.17,0.01))
  derivative <- diff(smooth) / diff(seq(0,7.17,0.01))
  y <- derivative
  x <- seq(0,7.17,0.01)[-1]
  z<-smooth[-1]
  
  data <- data.frame(x, y,z)
  data<-subset(data,!data$y=='NaN')
  rownames(data)<-1:nrow(data)
  
  baseline<-1.1*max(data$y[1:150])
  consecutive_up<-10
  
  # Find the point where the loess curve starts an upward trajectory
  change_point <- NULL
  for (i in 2:(nrow(data) - consecutive_up)) {
    upward_trend <- TRUE
    for (j in 1:consecutive_up) {
      if (!((data$y[i + j] - data$y[i + j - 1]) > 0.001)) {
        upward_trend <- FALSE
        break
      }
    }
    if (upward_trend && data$y[i] > baseline) {
      change_point <-data$x[i]
      break
    }
  }
  
  # Add a vertical line to indicate the point of upward trajectory
  if (!is.null(change_point) && tail(data$z, n = 1) > 1.1*head(data$z, n = 1)) {
    take_off.pt <- c(take_off.pt, change_point)
  } else {
    take_off.pt <- c(take_off.pt, NA)
  }
}

take_off.pt  # Return the vector of detected change points

take.off.pt<-data.frame(colnames(vsd_bgc)[1:41005],take_off.pt)
colnames(take.off.pt)<-c('Gene','Take-off point')

saveRDS(take.off.pt,"./take.off.pt.RDS")

######### plot all take off point ######

metadata_pst <- metadata_sub
metadata_pst$sampleID<-rownames(metadata_pst)
metadata_pst<-metadata_pst %>%
  select(sampleID,condition,PC1,PC2,UMAP1,UMAP2,pst_GMM.PCA)
metadata_pst$new.pst<-max(metadata_pst$pst_GMM.PCA)-metadata_pst$pst_GMM.PCA

take.off.pt<-subset(take.off.pt,!is.na(take.off.pt$`Take-off point`))
colnames(take.off.pt)<-c('Gene','new.pst')

gene_df<-as.data.frame(table(take.off.pt$new.pst))


TOP_df<-metadata_pst %>%
  dplyr::select(condition2,NMFrisk,new.pst)%>%
  subset(!is.na(new.pst))


density.plot<-ggplot(data=TOP_df,aes(x=new.pst,y=condition2,fill=condition2))+
  geom_density_ridges2(alpha=0.7)+
  theme_classic()+
  scale_fill_manual(values=color4disease)+
  coord_flip()

dev.new()
pdf('density.plot2.pdf')
density.plot
dev.off()

take.off.plot<-ggplot(data=take.off.pt,aes(x=new.pst))+
  geom_histogram(binwidth = 0.05,fill='black')+
  theme_classic()+
  scale_x_continuous(limits = c(0,9),breaks = seq(0,9,3))+
  geom_vline(xintercept=2)

dev.new()
pdf('take.off.plot2.pdf')
take.off.plot
dev.off()



