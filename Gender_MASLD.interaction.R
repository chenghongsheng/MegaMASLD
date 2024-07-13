setwd('./interaction')

library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(pheatmap)
library(limma)
library(ViSEAGO)
library(Boruta)
library(viridis)

# load dds data from 
dds<-readRDS("../Deseq2_analysis/dds.rds")

###########color setting###########
color4gender<-c("blue",'pink')

pie(rep(1,4),col=brewer.pal(8,"Set3")[c(7,1,6,4)])
color4disease<-c(brewer.pal(8,"Set3")[c(7,1,6)],"#FF6347")

###########DEGs###########
Hsa.dataset<-useDataset('hsapiens_gene_ensembl',mart=useMart("ensembl"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name',"gene_biotype"), 
               filters='ensembl_gene_id',
               values=rownames(dds),mart=Hsa.dataset)
genesymbols <- tapply(Genemap$external_gene_name, 
                      Genemap$ensembl_gene_id, paste, collapse="; ")

resultsNames(dds)

######## interaction of NAFL vs healthy ##############
res.NAFLvCtrl.MvF <- as.data.frame(results(dds, name="GenderM.conditionNAFL",alpha=0.05)) # condition effect (NAFLvshealthy) is different between genders or not
res.NAFLvCtrl.MvF<-res.NAFLvCtrl.MvF[order(res.NAFLvCtrl.MvF$padj),]
res.NAFLvCtrl.MvF$symbol<-genesymbols[rownames(res.NAFLvCtrl.MvF)]
sig.NAFLvCtrl.MvF<-res.NAFLvCtrl.MvF[which(res.NAFLvCtrl.MvF$padj<0.05 & 
                                             (abs(res.NAFLvCtrl.MvF$log2FoldChange)>1) &
                                             res.NAFLvCtrl.MvF$baseMean >100),]

write.csv(sig.NAFLvCtrl.MvF,'sig.NAFLvCtrl.MvF.csv')
write.csv(res.NAFLvCtrl.MvF,'res.NAFLvCtrl.MvF.csv')


######## interaction of NASH vs healthy ##############
res.NASHvCtrl.MvF <- as.data.frame(results(dds, name="GenderM.conditionNASH",alpha=0.05)) # condition effect (NASHvshealthy) is different between genders or not
res.NASHvCtrl.MvF<-res.NASHvCtrl.MvF[order(res.NASHvCtrl.MvF$padj),]
res.NASHvCtrl.MvF$symbol<-genesymbols[rownames(res.NASHvCtrl.MvF)]
sig.NASHvCtrl.MvF<-res.NASHvCtrl.MvF[which(res.NASHvCtrl.MvF$padj<0.05 & 
                                             (abs(res.NASHvCtrl.MvF$log2FoldChange)>1) &
                                             res.NASHvCtrl.MvF$baseMean >100),]

write.csv(sig.NASHvCtrl.MvF,'sig.NASHvCtrl.MvF.csv')
write.csv(res.NASHvCtrl.MvF,'res.NASHvCtrl.MvF.csv')


######## interaction of cirrhosis vs healthy ##############
res.cirrhosisvCtrl.MvF <- as.data.frame(results(dds, name="GenderM.conditioncirrhosis",alpha=0.05)) # condition effect (cirrhosisvshealthy) is different between genders or not
res.cirrhosisvCtrl.MvF<-res.cirrhosisvCtrl.MvF[order(res.cirrhosisvCtrl.MvF$padj),]
res.cirrhosisvCtrl.MvF$symbol<-genesymbols[rownames(res.cirrhosisvCtrl.MvF)]
sig.cirrhosisvCtrl.MvF<-res.cirrhosisvCtrl.MvF[which(res.cirrhosisvCtrl.MvF$padj<0.05 & 
                                                       (abs(res.cirrhosisvCtrl.MvF$log2FoldChange)>1) &
                                                       res.cirrhosisvCtrl.MvF$baseMean >100),]

write.csv(sig.cirrhosisvCtrl.MvF,'sig.cirrhosisvCtrl.MvF.csv')
write.csv(res.cirrhosisvCtrl.MvF,'res.cirrhosisvCtrl.MvF.csv')


######## interaction of NASH vs NAFL ##############
res.NASHvNAFL.MvF <- as.data.frame(results(dds, contrast=list("GenderM.conditionNASH","GenderM.conditionNAFL"), alpha=0.05))
res.NASHvNAFL.MvF<-res.NASHvNAFL.MvF[order(res.NASHvNAFL.MvF$padj),]
res.NASHvNAFL.MvF$symbol<-genesymbols[rownames(res.NASHvNAFL.MvF)]
sig.NASHvNAFL.MvF<-res.NASHvNAFL.MvF[which(res.NASHvNAFL.MvF$padj<0.05 & 
                                             (abs(res.NASHvNAFL.MvF$log2FoldChange)>1) &
                                             res.NASHvNAFL.MvF$baseMean >100),]

write.csv(sig.NASHvNAFL.MvF,'sig.NASHvNAFL.MvF.csv')
write.csv(res.NASHvNAFL.MvF,'res.NASHvNAFL.MvF.csv')


######## interaction of cirrhosis vs NAFL ##############
res.cirrhosisvNAFL.MvF <- as.data.frame(results(dds, contrast=list("GenderM.conditioncirrhosis","GenderM.conditionNAFL"), alpha=0.05))
res.cirrhosisvNAFL.MvF<-res.cirrhosisvNAFL.MvF[order(res.cirrhosisvNAFL.MvF$padj),]
res.cirrhosisvNAFL.MvF$symbol<-genesymbols[rownames(res.cirrhosisvNAFL.MvF)]
sig.cirrhosisvNAFL.MvF<-res.cirrhosisvNAFL.MvF[which(res.cirrhosisvNAFL.MvF$padj<0.05 & 
                                                       (abs(res.cirrhosisvNAFL.MvF$log2FoldChange)>1) &
                                                       res.cirrhosisvNAFL.MvF$baseMean >100),]

write.csv(sig.cirrhosisvNAFL.MvF,'sig.cirrhosisvNAFL.MvF.csv')
write.csv(res.cirrhosisvNAFL.MvF,'res.cirrhosisvNAFL.MvF.csv')

######## interaction of cirrhosis vs NASH ##############
res.cirrhosisvNASH.MvF <- as.data.frame(results(dds, contrast=list("GenderM.conditioncirrhosis","GenderM.conditionNASH"), alpha=0.05))
res.cirrhosisvNASH.MvF<-res.cirrhosisvNASH.MvF[order(res.cirrhosisvNASH.MvF$padj),]
res.cirrhosisvNASH.MvF$symbol<-genesymbols[rownames(res.cirrhosisvNASH.MvF)]
sig.cirrhosisvNASH.MvF<-res.cirrhosisvNASH.MvF[which(res.cirrhosisvNASH.MvF$padj<0.05 & 
                                                       (abs(res.cirrhosisvNASH.MvF$log2FoldChange)>1) &
                                                       res.cirrhosisvNASH.MvF$baseMean >100),]

write.csv(sig.cirrhosisvNASH.MvF,'sig.cirrhosisvNASH.MvF.csv')
write.csv(res.cirrhosisvNASH.MvF,'res.cirrhosisvNASH.MvF.csv')


#all DEGs
allDEG<-c(rownames(sig.NAFLvCtrl.MvF),rownames(sig.NASHvCtrl.MvF),rownames(sig.cirrhosisvCtrl.MvF),
          rownames(sig.NASHvNAFL.MvF),rownames(sig.cirrhosisvNAFL.MvF),rownames(sig.cirrhosisvNASH.MvF))
allDEG<-allDEG[!duplicated(allDEG)]


######## heatmap construction ##############
coldat_hm<-as.data.frame(dds@colData[,c(1,6)])
coldat_hm<-coldat_hm[order(coldat_hm$Gender,coldat_hm$condition),]

hm_mat<-assay(vsd_bc)
hm_mat<-subset(hm_mat,rownames(hm_mat) %in% allDEG)
hm_mat<-hm_mat[,rownames(coldat_hm)]

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

hm_mat<-hm_mat[,rownames(coldat_hm)]

break_hm = seq(-2.5, 2.5,length.out=100)

my_color_annotation<-list(condition= c("healthy"="#B3DE69","NAFL"="#8DD3C7",
                                       "NASH"="#FDB462","cirrhosis"="#FF6347"),
                          Gender=c("M"="blue",'F'='pink'))


pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         angle_col = 45,gaps_col = 411,cutree_row=2)



hm_mat2<-assay(vsd_bgc)
hm_mat2<-subset(hm_mat2,rownames(hm_mat2) %in% allDEG)
hm_mat2<-hm_mat2[,rownames(coldat_hm)]

pheatmap(hm_mat2,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         angle_col = 45,gaps_col = 411)

######## consistently up or down DEGs ##############
sig.NAFLvCtrl.MvF.up<-subset(sig.NAFLvCtrl.MvF,sig.NAFLvCtrl.MvF$log2FoldChange>0)
sig.NASHvCtrl.MvF.up<-subset(sig.NASHvCtrl.MvF,sig.NASHvCtrl.MvF$log2FoldChange>0)
sig.cirrhosisvCtrl.MvF.up<-subset(sig.cirrhosisvCtrl.MvF,sig.cirrhosisvCtrl.MvF$log2FoldChange>0)
sig.NASHvNAFL.MvF.up<-subset(sig.NASHvNAFL.MvF,sig.NASHvNAFL.MvF$log2FoldChange>0)
sig.cirrhosisvNAFL.MvF.up<-subset(sig.cirrhosisvNAFL.MvF,sig.cirrhosisvNAFL.MvF$log2FoldChange>0)
sig.cirrhosisvNASH.MvF.up<-subset(sig.cirrhosisvNASH.MvF,sig.cirrhosisvNASH.MvF$log2FoldChange>0)

allDEG.up<-intersect(rownames(sig.NAFLvCtrl.MvF.up),rownames(sig.NASHvCtrl.MvF.up))
allDEG.up<-intersect(allDEG.up,rownames(sig.cirrhosisvCtrl.MvF.up))

#allDEG.up<-allDEG.up[!duplicated(allDEG.up)]


sig.NAFLvCtrl.MvF.down<-subset(sig.NAFLvCtrl.MvF,sig.NAFLvCtrl.MvF$log2FoldChange<0)
sig.NASHvCtrl.MvF.down<-subset(sig.NASHvCtrl.MvF,sig.NASHvCtrl.MvF$log2FoldChange<0)
sig.cirrhosisvCtrl.MvF.down<-subset(sig.cirrhosisvCtrl.MvF,sig.cirrhosisvCtrl.MvF$log2FoldChange<0)
sig.NASHvNAFL.MvF.down<-subset(sig.NASHvNAFL.MvF,sig.NASHvNAFL.MvF$log2FoldChange<0)
sig.cirrhosisvNAFL.MvF.down<-subset(sig.cirrhosisvNAFL.MvF,sig.cirrhosisvNAFL.MvF$log2FoldChange<0)
sig.cirrhosisvNASH.MvF.down<-subset(sig.cirrhosisvNASH.MvF,sig.cirrhosisvNASH.MvF$log2FoldChange<0)

allDEG.down<-intersect(rownames(sig.NAFLvCtrl.MvF.down),rownames(sig.NASHvCtrl.MvF.down))
allDEG.down<-intersect(allDEG.down,rownames(sig.cirrhosisvCtrl.MvF.down))

intersect(allDEG.up,allDEG.down) #NULL


sig.NAFLvCtrl.MvF.sel.up<-subset(sig.NAFLvCtrl.MvF,rownames(sig.NAFLvCtrl.MvF) %in% allDEG.up)
sig.NASHvCtrl.MvF.sel.up<-subset(sig.NASHvCtrl.MvF,rownames(sig.NASHvCtrl.MvF) %in% allDEG.up)
sig.cirrhosisvCtrl.MvF.sel.up<-subset(sig.cirrhosisvCtrl.MvF, rownames(sig.cirrhosisvCtrl.MvF) %in% allDEG.up)


sig.NAFLvCtrl.MvF.sel.down<-subset(sig.NAFLvCtrl.MvF,rownames(sig.NAFLvCtrl.MvF) %in% allDEG.down)
sig.NASHvCtrl.MvF.sel.down<-subset(sig.NASHvCtrl.MvF,rownames(sig.NASHvCtrl.MvF) %in% allDEG.down)
sig.cirrhosisvCtrl.MvF.sel.down<-subset(sig.cirrhosisvCtrl.MvF, rownames(sig.cirrhosisvCtrl.MvF) %in% allDEG.down)

write.csv(sig.NAFLvCtrl.MvF.sel.up,'sig.NAFLvCtrl.MvF.sel.up.csv')
write.csv(sig.NASHvCtrl.MvF.sel.up,'sig.NASHvCtrl.MvF.sel.up.csv')
write.csv(sig.cirrhosisvCtrl.MvF.sel.up,'sig.cirrhosisvCtrl.MvF.sel.up.csv')

write.csv(sig.NAFLvCtrl.MvF.sel.down,'sig.NAFLvCtrl.MvF.sel.down.csv')
write.csv(sig.NASHvCtrl.MvF.sel.down,'sig.NASHvCtrl.MvF.sel.down.csv')
write.csv(sig.cirrhosisvCtrl.MvF.sel.down,'sig.cirrhosisvCtrl.MvF.sel.down.csv')

#heatmap
allDEG.updown<-c(allDEG.up,allDEG.down)

hm_mat<-subset(hm_mat,rownames(hm_mat) %in% allDEG.updown)

pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         angle_col = 45,gaps_col = 411,cutree_row=2)

hm_mat2<-subset(hm_mat2,rownames(hm_mat2) %in% allDEG.updown)

pheatmap(hm_mat2,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         angle_col = 45,gaps_col = 411,cutree_row=3)


allDEG.up<-as.data.frame(allDEG.up)
allDEG.up$symbols<-genesymbols[allDEG.up$allDEG.up]

allDEG.down<-as.data.frame(allDEG.down)
allDEG.down$symbols<-genesymbols[allDEG.down$allDEG.down]

write.csv(allDEG.up,'allDEG.up.csv')
write.csv(allDEG.down,'allDEG.down.csv')


######## DEGs Gender_all samples  ##############
res.MvF <- as.data.frame(results(dds, contrast=c('Gender','M','F'),alpha=0.05)) # gender effect
res.MvF<-res.MvF[order(res.MvF$padj),]
res.MvF$symbol<-genesymbols[rownames(res.MvF)]
sig.MvF<-res.MvF[which(res.MvF$padj<0.05 &
                         (abs(res.MvF$log2FoldChange)>1) &
                         res.MvF$baseMean >100),] 

#volcanoplot
res.MvF <-res.MvF %>%
  mutate(intgene=case_when(rownames(res.MvF) %in% allDEG.up$allDEG.up ~ 'up',
                           rownames(res.MvF) %in% allDEG.down$allDEG.down ~ 'down',
                           T~"not"))


res.MvF <-res.MvF %>%
  mutate(lab=case_when(intgene %in% c('up','down') ~symbol,
                       T~""))

res.MvF$intgene<-factor(res.MvF$intgene,levels = c('down','up','not'))
res.MvF<-res.MvF[rev(order(res.MvF$intgene)),]

intgene.volcano<-ggplot(res.MvF,aes(x=log2FoldChange,y=-log10(padj),color=intgene,alpha=intgene))+
  geom_point(size=2)+
  scale_alpha_manual(values=c(1,1,0.2))+
  scale_color_manual(values=c('red','red','grey70'))+
  geom_text_repel(aes(label=lab),max.overlaps = 50,color='black')+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(0,75))+
  theme_light()

dev.new()
pdf("intgene.volcano.pdf")
intgene.volcano
dev.off()

# violin plot of selected genes
norm_count<-counts(dds,normalized=T)

sel.norm_count<-as.data.frame(t(subset(norm_count,rownames(norm_count) %in% c("ENSG00000012817",'ENSG00000114374','ENSG00000198692'))))

colnames(sel.norm_count)<-genesymbols[c("ENSG00000012817",'ENSG00000114374','ENSG00000198692')]
sel.norm_count<-sel.norm_count[rownames(coldat_hm),]
all(rownames(sel.norm_count)==rownames(coldat_hm))

sel.norm_count$Gender<-coldat_hm$Gender
sel.norm_count$condition<-coldat_hm$condition

KDM5D<-ggplot(sel.norm_count,aes(x=Gender,y=KDM5D,fill=condition))+
  geom_violin()+
  scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values=color4disease)+
  stat_summary(fun = 'mean',geom='crossbar',position = "dodge")+
  ggtitle("bar is mean; KDM5D")+
  theme_light()

dev.new()
pdf("KDM5D.pdf")
KDM5D
dev.off()


USP9Y<-ggplot(sel.norm_count,aes(x=Gender,y=USP9Y,fill=condition))+
  geom_violin()+
  scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values=color4disease)+
  stat_summary(fun = 'mean',geom='crossbar',position = "dodge")+
  ggtitle("bar is mean; USP9Y")+
  theme_light()

dev.new()
pdf("USP9Y.pdf")
USP9Y
dev.off()


EIF1AY<-ggplot(sel.norm_count,aes(x=Gender,y=EIF1AY,fill=condition))+
  geom_violin()+
  scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values=color4disease)+
  stat_summary(fun = 'mean',geom='crossbar',position = "dodge")+
  ggtitle("bar is mean; EIF1AY")+
  theme_light()

dev.new()
pdf("EIF1AY.pdf")
EIF1AY
dev.off()
