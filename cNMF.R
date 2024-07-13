setwd('./cNMF')

library(NMF)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(limma)
library(viridis)
library(ViSEAGO)

######## load data ##############
color4gender<-c("blue",'pink')

pie(rep(1,4),col=brewer.pal(8,"Set3")[c(7,1,6,4)])
color4disease<-c(brewer.pal(8,"Set3")[c(7,1,6)],"#FF6347")

metadata<-read.csv("../Deseq2_analysis/metadata.csv",header=T,check.names = F,row.names = 1)
metadata$condition<- factor(metadata$condition,levels=c("healthy",'NAFL','NASH','cirrhosis'))
metadata$batch<-factor(metadata$batch)
metadata$pred.gender<-factor(metadata$pred.gender,levels=c("M","F"))

vsd.bgc<-read.table("../Deseq2_analysis/vsd.bgc_ssGSEA.txt",header=T,check.names = F) # bg corrected

background<-read.csv("../Deseq2_analysis/expressed_genes.txt", sep="")

#vsd of expressed genes
vsd.bgc.sel<-subset(vsd.bgc,rownames(vsd.bgc) %in% background$x)
vsd.bgc.sel<-vsd.bgc.sel[!duplicated(vsd.bgc.sel$symbols),] #remove duplicated symbols
vsd.bgc.sel<-subset(vsd.bgc.sel,!(vsd.bgc.sel$symbols==""|is.na(vsd.bgc.sel$symbols))) #remove empty gene symbol or NA

rownames(vsd.bgc.sel)<-vsd.bgc.sel$symbols
vsd.bgc.sel$symbols<-NULL

vsd.bgc.sel$var<-rowVars(as.matrix(vsd.bgc.sel),useNames = F)
vsd.bgc.sel<-vsd.bgc.sel[order(vsd.bgc.sel$var,decreasing = T),]
vsd.bgc.top1000<-vsd.bgc.sel[1:1000,] #select top1000 more variable
vsd.bgc.top1000$var<-NULL
vsd.bgc.top1000<-vsd.bgc.top1000+abs(min(vsd.bgc.top1000))+0.1

#rm(vsd.bgc)
#rm(background)

######## NMF ##############
set.seed(123456)

#rank 1-49
res <- nmf(vsd.bgc.top1000, seq(1,50,2), method='brunet', nrun=10, seed=123456)
plot(res)
consensusmap(res$consensus[2:10])


dev.new()
pdf("NMF_rank1-49.pdf",width=15)
plot(res)
dev.off()

dev.new()
pdf("NMF_consensusmap_rank3.pdf")
consensusmap(res$consensus[2])
dev.off() #narrow to rank 2-5

dev.new()
pdf("NMF_consensusmap_rank9.pdf")
consensusmap(res$consensus[5])
dev.off() 

saveRDS(res,"NMF_res1-49.RDS")

#rank 2-5 
res2 <- nmf(vsd.bgc.top1000, seq(2,5), method='brunet', nrun=10, seed=123456)
plot(res2)
consensusmap(res2)


dev.new()
pdf("NMF_rank2-5.pdf",width=15)
plot(res2)
dev.off()

dev.new()
pdf("NMF_consensusmap_rank2-5.pdf",width=10,height = 10)
consensusmap(res2)
dev.off() 


#use rank 2
res_rank2 <- nmf(vsd.bgc.top1000, 2, method='brunet', nrun=10, seed=123456)
plot(res_rank2)
consensusmap(res_rank2)


NMF.feature<-as.data.frame(basis(res_rank2))
NMF.sample<-coef(res_rank2)
s<-featureScore(res_rank2) 

summary(s)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.01070 0.03901 0.05649 0.08039 0.50480 

s<-extractFeatures(res_rank2) 
str(s)

# obtain top 50 differential features
colnames(NMF.feature)<-c('gs1','gs2')

NMF.feature$dif<-NMF.feature$gs1-NMF.feature$gs2
NMF.feature<-NMF.feature[order(NMF.feature$dif,decreasing = T),]
NMF.feature<-NMF.feature %>%
  mutate(clusfeature=case_when(dif>0 ~ 'gs1',
                               dif<0 ~ "gs2"))

clus1_feature<-NMF.feature[1:50,]
clus2_feature<-NMF.feature[951:1000,]



#plot(reorder(1:length(NMF.feature$dif),NMF.feature$dif,decreasing = T),NMF.feature$dif)


# obtain clustering results of samples
NMF.sample<-as.data.frame(t(NMF.sample))
colnames(NMF.sample)<-c('score1','score2')
sil<-as.data.frame(silhouette(res_rank2))

all(rownames(sil)==rownames(NMF.sample))
NMF.sample$cluster<-sil$cluster
NMF.sample$dif<-NMF.sample$score1-NMF.sample$score2

#add cluster to metadata
all(rownames(metadata)==rownames(NMF.sample))
metadata$NMFclus<-NMF.sample$cluster
metadata$NMFclus<-factor(metadata$NMFclus)
metadata$NMFdiff<-NMF.sample$dif


######## heatmap based on NMF ##############
hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
break_hm = seq(-2, 2,length.out=100)

my_color_annotation<-list(condition= c("healthy"="#CFF196","NAFL"="#8DD3C7",
                                       "NASH"="#FDB462","cirrhosis"="#FF6347"),
                          NMFclus=c('1'='red','2'='black'))

vsd_mat_NMF100 <-vsd.bgc.top1000[c(rownames(clus1_feature),rownames(clus2_feature)),rownames(NMF.sample[order(NMF.sample$dif,decreasing = F),])]

heatmap.top100.NMF<-pheatmap(vsd_mat_NMF100,scale="row",border_color = NA,color = hm_color,
                             show_rownames = F,show_colnames = F,
                             cluster_rows = T,cluster_cols =F,
                             annotation_col = coldat_hm,
                             breaks = break_hm,
                             annotation_colors = my_color_annotation,
                             clustering_distance_rows = "correlation",
                             clustering_distance_cols = "correlation",
                             angle_col = 45) 

dev.new()
pdf("heatmap.top100.NMF.pdf")
heatmap.top100.NMF
dev.off() 


#export data from NMF + silhouette
write.csv(metadata,"metadata_NMF.csv")
write.csv(NMF.feature[c(rownames(clus1_feature),rownames(clus2_feature)),],"NMF.markers.csv")

######## viseago of NMF gene set ##############
clus1_ensembl<-rownames(subset(vsd.bgc,vsd.bgc$symbols %in% rownames(clus1_feature)))
clus2_ensembl<-rownames(subset(vsd.bgc,vsd.bgc$symbols %in% rownames(clus2_feature)))

background<-scan("../Deseq2_analysis/expressed_genes.txt",
                 quiet=TRUE,
                 what="")

Ensembl<-ViSEAGO::Ensembl2GO()
ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "hsapiens_gene_ensembl",
  Ensembl)

#clus1 gene
BP_clus1_genes<-ViSEAGO::create_topGOdata(
  geneSel=clus1_ensembl,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus1_genes<-topGO::runTest(
  BP_clus1_genes,
  algorithm ="classic",
  statistic = "fisher")

BP_clus1_genes_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_clus1_genes","classic_clus1_genes")))

BP_clus1_genes_res<-as.data.frame(BP_clus1_genes_sResults@data)
write.table(BP_clus1_genes_res,"BP_clus1_genes.NMF_res.txt")


#clus2 gene
BP_clus2_genes<-ViSEAGO::create_topGOdata(
  geneSel=clus2_ensembl,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus2_genes<-topGO::runTest(
  BP_clus2_genes,
  algorithm ="classic",
  statistic = "fisher")

BP_clus2_genes_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_clus2_genes","classic_clus2_genes")))

BP_clus2_genes_res<-as.data.frame(BP_clus2_genes_sResults@data)
write.table(BP_clus2_genes_res,"BP_clus2_genes.NMF_res.txt")


###########NMF of MASL livers ###########
metadata <- metadata %>%
  mutate(NMFrisk=case_when(NMFclus=="1"~"high",
                           NMFclus=="2"~"low"))
metadata$NMFrisk<-factor(metadata$NMFrisk,levels=c("low",'high'))

###########import dds data ###########

dds<-readRDS("../Deseq2_analysis/dds.rds")

count<-dds@assays@data$counts
count<-as.data.frame(count)
#rm(dds)

all(colnames(count)==rownames(metadata)) #TRUE


###########subset NAFL ###########
metadata.NAFL<-subset(metadata,metadata$condition=="NAFL")
metadata.NAFL$condition<-as.character(metadata.NAFL$condition)
coldat.NAFL<-metadata.NAFL %>% dplyr::select("condition",'batch','pred.gender','NMFrisk')

count.NAFL<-count[,rownames(metadata.NAFL)]

dds.NAFL<-DESeqDataSetFromMatrix(countData = count.NAFL,colData =coldat.NAFL,design =~batch+pred.gender+NMFrisk)
dds.NAFL<-DESeq(dds.NAFL)

vsd.NAFL<-varianceStabilizingTransformation(dds.NAFL, blind=TRUE)

vsd.NAFL_bgc<-vsd.NAFL
vsd.NAFL_bgcmat<-assay(vsd.NAFL_bgc)
mm<-model.matrix(~NMFrisk,colData(vsd.NAFL))
vsd.NAFL_bgcmat<-removeBatchEffect(vsd.NAFL_bgcmat,batch=vsd.NAFL_bgc$batch,batch2 = vsd.NAFL_bgc$pred.gender,design=mm)
assay(vsd.NAFL_bgc)<-vsd.NAFL_bgcmat


pcadat_vsd.NAFL_bgc<-plotPCA(vsd.NAFL_bgc,intgroup="NMFrisk",returnData=T)
plotPCA(vsd.NAFL_bgc,intgroup="pred.gender")

all(rownames(pcadat_vsd.NAFL_bgc)==rownames(coldat)) #sanity check

pcadat_vsd.NAFL_bgc$batch<-coldat.NAFL$batch
pcadat_vsd.NAFL_bgc$gender<-coldat.NAFL$pred.gender

percentVar.vsd.NAFL.bgc<-round(100*attr(pcadat_vsd.NAFL_bgc,"percentVar"))

ggplot(pcadat_vsd.NAFL_bgc, aes(PC1, PC2, fill=gender)) +
  geom_point(size=3,pch=21,stroke=0.5,alpha=0.9)+ 
  #geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd.NAFL.bgc[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd.NAFL.bgc[2],"% variance")) + 
  #stat_ellipse()+
  # geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = color4gender)+
  scale_color_manual(values = color4gender)+
  coord_fixed()


pca.NMFrisk_bgc<-ggplot(pcadat_vsd.NAFL_bgc, aes(PC1, PC2, fill=NMFrisk)) +
  geom_point(size=3,pch=21,stroke=0.5,alpha=0.9)+ 
  #geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd.NAFL.bgc[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd.NAFL.bgc[2],"% variance")) + 
  #stat_ellipse()+
  # geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("black",'red'))+
  scale_color_manual(values =c("black",'red'))+
  coord_fixed()

dev.new()
pdf("pca.NMFrisk_bgc.pdf")
pca.NMFrisk_bgc
dev.off()

###########DEGs###########
Hsa.dataset<-useDataset('hsapiens_gene_ensembl',mart=useMart("ensembl"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name',"gene_biotype"), 
               filters='ensembl_gene_id',
               values=rownames(dds.NAFL),mart=Hsa.dataset)
genesymbols <- tapply(Genemap$external_gene_name, 
                      Genemap$ensembl_gene_id, paste, collapse="; ")

######## NMFrisk ##############
######## High vs Low ##############
res.highvlow<-as.data.frame(results(dds.NAFL,contrast=c('NMFrisk','high','low'),alpha=0.05))
res.highvlow<-res.highvlow[order(res.highvlow$padj),]
res.highvlow$symbol<-genesymbols[rownames(res.highvlow)]
sig.highvlow<-res.highvlow[which((res.highvlow$pvalue<0.05)& (abs(res.highvlow$log2FoldChange)>1) & res.highvlow$baseMean >100),]
write.csv(sig.highvlow,"sig.highvlow_bgc.csv")
write.csv(res.highvlow,"res.highvlow_bgc.csv")

########### demographic - gender ###########
df_gender<-metadata.NAFL %>%
  group_by(pred.gender,NMFrisk) %>%
  summarise(total=length(NMFrisk))

write.csv(df_gender,"df_gender.forchisquare.csv")

M_pie<-ggplot(subset(df_gender,df_gender$pred.gender=="M"),aes(x="",y=total,fill=NMFrisk))+
  geom_bar(stat = "identity",width=1,color="white")+
  scale_fill_manual(values = c("black",'red'))+
  coord_polar("y",start=0)+
  theme_void()+
  geom_text(aes(y = c(180,60), label = total), color = "white", size=6)

F_pie<-ggplot(subset(df_gender,df_gender$pred.gender=="F"),aes(x="",y=total,fill=NMFrisk))+
  geom_bar(stat = "identity",width=1,color="white")+
  scale_fill_manual(values = c("black",'red'))+
  coord_polar("y",start=0)+
  theme_void()+
  geom_text(aes(y = c(160,30), label = total), color = "white", size=6)


dev.new()
pdf("M_pie.pdf")
M_pie
dev.off()


dev.new()
pdf("F_pie.pdf")
F_pie
dev.off()

########### demographic - age ###########
t.test(Age~NMFrisk,data=metadata.NAFL)

df_age<- subset(metadata.NAFL,!is.na(metadata.NAFL$Age)) %>%
  group_by(NMFrisk) %>%
  summarise(mean.age=mean(Age),sd.age=sd(Age),median.age=median(Age))
write.csv(df_age,'df_age.csv') #need descriptive data

AgeNMF<-ggplot(metadata.NAFL,aes(x=NMFrisk,y=Age,fill=NMFrisk))+
  geom_violin(alpha=0.2)+
  scale_fill_manual(values = c("black",'red'))+
  geom_dotplot(stackdir = "center",binaxis = "y",binwidth = 1.5)+
  stat_summary(fun="mean",geom="crossbar",width=0.8)+
  #stat_summary(fun="median",geom="crossbar",width=0.8,color="grey80")
  ggtitle("bar is mean; lo=50.42, hi=57.09; p= 6.889e-05")

dev.new()
pdf("AgeNMF.pdf")
AgeNMF
dev.off()

df_age.bygender<- subset(metadata.NAFL,!is.na(metadata.NAFL$Age)) %>%
  group_by(pred.gender,NMFrisk) %>%
  summarise(mean.age=mean(Age))

metadata.NAFL%>%
  #  subset(metadata.NAFL,!is.na(metadata.NAFL$Age)) %>%
  group_by(pred.gender) %>%
  pairwise_t_test(Age ~ NMFrisk,p.adjust.method = "holm")

AgeNMF.bygender<-ggplot(metadata.NAFL,aes(x=NMFrisk,y=Age,fill=NMFrisk))+
  geom_violin(alpha=0.2)+
  scale_fill_manual(values = c("black",'red'))+
  geom_dotplot(stackdir = "center",binaxis = "y",binwidth = 1.5)+
  facet_grid(~pred.gender)+
  stat_summary(fun="mean",geom="crossbar",width=0.8)+
  #stat_summary(fun="median",geom="crossbar",width=0.8,color="grey80")
  ggtitle("bar is mean; lo.M=49.6, hi.M=56.2;lo.F=51.5, hi.F=58.7; p.M= 0.00213, p.F=0.00614")

dev.new()
pdf("AgeNMF.bygender.pdf",width=10)
AgeNMF.bygender
dev.off()


########### demographic - BMI ########### 
t.test(BMI~NMFrisk,data=metadata.NAFL)# not enough data


########### selected genes ########### 
vsd.NAFL_bgc.df<-as.data.frame(vsd.NAFL_bgcmat)

vsd.NAFL_bgc.df$symbol<-genesymbols[rownames(vsd.NAFL_bgc.df)]
vsd.NAFL_bgc.df<-subset(vsd.NAFL_bgc.df,!duplicated(vsd.NAFL_bgc.df$symbol))
vsd.NAFL_bgc.df <-subset(vsd.NAFL_bgc.df,!vsd.NAFL_bgc.df$symbol=="")

rownames(vsd.NAFL_bgc.df)<-vsd.NAFL_bgc.df$symbol
vsd.NAFL_bgc.df$symbol<-NULL

all(colnames(vsd.NAFL_bgc.df)==rownames(metadata.NAFL)) #sanity check

vsd.NAFL_bgc.df<-as.data.frame(t(vsd.NAFL_bgc.df))
all(rownames(vsd.NAFL_bgc.df)==rownames(metadata.NAFL)) #sanity check
vsd.NAFL_bgc.df$NMFrisk<-metadata.NAFL$NMFrisk


gene='IGHA1'

IGHA1.plot<-ggplot(vsd.NAFL_bgc.df,aes(x=NMFrisk,y=IGHA1,fill=NMFrisk))+
  geom_violin()+
  scale_fill_manual(values=c('black','red'))+
  stat_summary(fun=mean,geom="crossbar",color="white",width=0.5)+
  labs(fill="",x="",y='Normalized expression')+
  ggtitle(paste(as.character(gene),' by NMF stratification'))+
  theme_light()+
  theme(axis.text.x.bottom=element_text(size = rel(2),colour = 'black'),
        axis.text.y.left = element_text(size=rel(2),color = 'black'),
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(2),face = 'bold'),
        legend.position = "none")

dev.new()
pdf("IGHA1.MASL.pdf",width=10)
IGHA1.plot
dev.off()
