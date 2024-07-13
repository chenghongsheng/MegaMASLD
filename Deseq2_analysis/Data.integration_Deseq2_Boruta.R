setwd('./Deseq2_analysis')

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

###########color setting###########
color4gender<-c("blue",'pink')

pie(rep(1,4),col=brewer.pal(8,"Set3")[c(7,1,6,4)])
color4disease<-c(brewer.pal(8,"Set3")[c(7,1,6)],"#FF6347")

###############DESeq2 batch correction################

dds<-DESeqDataSetFromMatrix(countData=Fcount_all,colData = coldata_all, design=~batch + Gender + condition)
dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

saveRDS(dds,'dds.RDS')

#PCA plot without batch and gender correction
plotPCA(vsd,intgroup="condition") 

#PCA plot after batch correction
vsd_bc<-vsd
vsd_bcmat<-assay(vsd_bc)
mm<-model.matrix(~Gender+condition,colData(vsd))
vsd_bcmat<-removeBatchEffect(vsd_bcmat,batch=vsd_bc$batch,design=mm)
assay(vsd_bc)<-vsd_bcmat

plotPCA(vsd_bc,intgroup="Gender")

#PCA plot after batch and gender correction
vsd_bgc<-vsd
vsd_bgcmat<-assay(vsd_bgc)
mm<-model.matrix(~condition,colData(vsd))
vsd_bgcmat<-removeBatchEffect(vsd_bgcmat,batch=vsd_bgc$batch,batch2 = vsd_bgc$Gender,design=mm)
assay(vsd_bgc)<-vsd_bgcmat

plotPCA(vsd_bgc,intgroup="condition")
plotPCA(vsd_bgc,intgroup="Gender")



########### DEGs ###########
Hsa.dataset<-useDataset('hsapiens_gene_ensembl',mart=useMart("ensembl"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name',"gene_biotype"), 
               filters='ensembl_gene_id',
               values=rownames(dds),mart=Hsa.dataset)
genesymbols <- tapply(Genemap$external_gene_name, 
                      Genemap$ensembl_gene_id, paste, collapse="; ")

###### NAFL vs healthy ##############
res.NAFLvCtrl<-as.data.frame(results(dds,contrast=c('condition','NAFL','healthy'),alpha=0.05))
res.NAFLvCtrl<-res.NAFLvCtrl[order(res.NAFLvCtrl$padj),]
res.NAFLvCtrl$symbol<-genesymbols[rownames(res.NAFLvCtrl)]
sig.NAFLvCtrl<-res.NAFLvCtrl[which((res.NAFLvCtrl$padj<0.05)& (abs(res.NAFLvCtrl$log2FoldChange)>1) & res.NAFLvCtrl$baseMean >100),]
write.csv(sig.NAFLvCtrl,"sig.NAFLvCtrl_bgc.csv")
write.csv(res.NAFLvCtrl,"res.NAFLvCtrl_bgc.csv")


###### NASH vs healthy ##############
res.NASHvCtrl<-as.data.frame(results(dds,contrast=c('condition','NASH','healthy'),alpha=0.05))
res.NASHvCtrl<-res.NASHvCtrl[order(res.NASHvCtrl$padj),]
res.NASHvCtrl$symbol<-genesymbols[rownames(res.NASHvCtrl)]
sig.NASHvCtrl<-res.NASHvCtrl[which((res.NASHvCtrl$padj<0.05)& (abs(res.NASHvCtrl$log2FoldChange)>1) & res.NASHvCtrl$baseMean >100),]
write.csv(sig.NASHvCtrl,"sig.NASHvCtrl_bgc.csv")
write.csv(res.NASHvCtrl,"res.NASHvCtrl_bgc.csv")

###### cirrhosis vs healthy ##############
res.cirrhosisvCtrl<-as.data.frame(results(dds,contrast=c('condition','cirrhosis','healthy'),alpha=0.05))
res.cirrhosisvCtrl<-res.cirrhosisvCtrl[order(res.cirrhosisvCtrl$padj),]
res.cirrhosisvCtrl$symbol<-genesymbols[rownames(res.cirrhosisvCtrl)]
sig.cirrhosisvCtrl<-res.cirrhosisvCtrl[which((res.cirrhosisvCtrl$padj<0.05)& (abs(res.cirrhosisvCtrl$log2FoldChange)>1) &
                                               res.cirrhosisvCtrl$baseMean >100),]
write.csv(sig.cirrhosisvCtrl,"sig.cirrhosisvCtrl_bgc.csv")
write.csv(res.cirrhosisvCtrl,"res.cirrhosisvCtrl_bgc.csv")

###### NASH vs NAFL ##############
res.NASHvNAFL<-as.data.frame(results(dds,contrast=c('condition','NASH','NAFL'),alpha=0.05))
res.NASHvNAFL<-res.NASHvNAFL[order(res.NASHvNAFL$padj),]
res.NASHvNAFL$symbol<-genesymbols[rownames(res.NASHvNAFL)]
sig.NASHvNAFL<-res.NASHvNAFL[which((res.NASHvNAFL$padj<0.05)& (abs(res.NASHvNAFL$log2FoldChange)>1) & res.NASHvNAFL$baseMean>100),]
write.csv(sig.NASHvNAFL,"sig.NASHvNAFL_bgc.csv")
write.csv(res.NASHvNAFL,"res.NASHvNAFL_bgc.csv")

###### cirrhosis vs NAFL ##############
res.cirrhosisvNAFL<-as.data.frame(results(dds,contrast=c('condition','cirrhosis','NAFL'),alpha=0.05))
res.cirrhosisvNAFL<-res.cirrhosisvNAFL[order(res.cirrhosisvNAFL$padj),]
res.cirrhosisvNAFL$symbol<-genesymbols[rownames(res.cirrhosisvNAFL)]
sig.cirrhosisvNAFL<-res.cirrhosisvNAFL[which((res.cirrhosisvNAFL$padj<0.05)& (abs(res.cirrhosisvNAFL$log2FoldChange)>1) &
                                               res.cirrhosisvNAFL$baseMean>100),]
write.csv(sig.cirrhosisvNAFL,"sig.cirrhosisvNAFL_bgc.csv")
write.csv(res.cirrhosisvNAFL,"res.cirrhosisvNAFL_bgc.csv")

###### cirrhosis vs NASH ##############
res.cirrhosisvNASH<-as.data.frame(results(dds,contrast=c('condition','cirrhosis','NASH'),alpha=0.05))
res.cirrhosisvNASH<-res.cirrhosisvNASH[order(res.cirrhosisvNASH$padj),]
res.cirrhosisvNASH$symbol<-genesymbols[rownames(res.cirrhosisvNASH)]
sig.cirrhosisvNASH<-res.cirrhosisvNASH[which((res.cirrhosisvNASH$padj<0.05)& (abs(res.cirrhosisvNASH$log2FoldChange)>1) & 
                                               res.cirrhosisvNASH$baseMean>100),]
write.csv(sig.cirrhosisvNASH,"sig.cirrhosisvNASH_bgc.csv")
write.csv(res.cirrhosisvNASH,"res.cirrhosisvNASH_bgc.csv")


###### plot.heatmap ##############
conditionDEG.all<-c(rownames(sig.NAFLvCtrl),rownames(sig.NASHvCtrl),rownames(sig.cirrhosisvCtrl),
                    rownames(sig.NASHvNAFL),rownames(sig.cirrhosisvNAFL),rownames(sig.cirrhosisvNASH))

conditionDEG.all<-conditionDEG.all[!duplicated(conditionDEG.all)]

hm_mat.bgc<-assay(vsd_bgc)
hm_condition<-subset(hm_mat.bgc,rownames(hm_mat.bgc) %in% conditionDEG.all)

coldat_hm<-as.data.frame(dds@colData[,c(1,6)])
coldat_hm<-coldat_hm[order(coldat_hm$condition),]

hm_condition<-hm_condition[,rownames(coldat_hm)]
hm_condition<-as.data.frame(hm_condition)

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
break_hm = seq(-2, 2,length.out=100)

my_color_annotation<-list(condition= c("healthy"="#B3DE69","NAFL"="#8DD3C7",
                                       "NASH"="#FDB462","cirrhosis"="#FF6347"))

hm.plot<-pheatmap(hm_condition,scale="row",border_color = NA,color = hm_color,
                  show_rownames = F,show_colnames = F,
                  cluster_rows = T,cluster_cols =F,
                  annotation_col = coldat_hm[,1,drop=F],
                  breaks = break_hm,
                  annotation_colors = my_color_annotation,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "correlation",
                  angle_col = 45) 

dev.new()
pdf("hm.plot.pdf")
hm.plot
dev.off()

####### Functional enrichment analysis ####

hc<-hclust(as.dist(1-cor(t(hm_condition),method = "pearson")),method = "complete")
all(hc$order==hm.plot$tree_row$order) #all TRUE
hm_cluster <- as.data.frame(cutree(tree = hc, k = 2))
table(hm_cluster)

NASH_genes<-subset(hm_cluster,hm_cluster$`cutree(tree = hc, k = 2)`=="1")
healthy_genes<-subset(hm_cluster,hm_cluster$`cutree(tree = hc, k = 2)`=="2")

expressed_genes<-as.data.frame(count)
expressed_genes[expressed_genes=="0"]<-NA
expressed_genes$NAcount<-rowSums(is.na(expressed_genes))
expressed_genes<-expressed_genes %>%
  mutate(selected=case_when(NAcount>415 ~ 'N',
                            T ~ 'Y'))
expressed_genes<-subset(expressed_genes,expressed_genes$selected=="Y")
expressed_genes<-rownames(expressed_genes)
write.table(expressed_genes,'expressed_genes.txt')

background<-scan("expressed_genes.txt",
                 quiet=TRUE,
                 what="")

Ensembl<-ViSEAGO::Ensembl2GO()
ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "hsapiens_gene_ensembl",
  Ensembl)


#NASH gene
write.table(rownames(NASH_genes),'NASHgene.txt')

NASH_genes<-scan(
  "NASHgene.txt",
  quiet=TRUE,
  what="")

BP_NASH_genes<<-ViSEAGO::create_topGOdata(
  geneSel=NASH_genes,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_NASH_genes<-topGO::runTest(
  BP_NASH_genes,
  algorithm ="classic",
  statistic = "fisher")

BP_NASH_genes_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_NASH_genes","classic_NASH_genes")))

BP_NASH_genes_res<-as.data.frame(BP_NASH_genes_sResults@data)
write.table(BP_NASH_genes_res,"BP_NASH_genes_res.txt")

#healthy gene
write.table(rownames(healthy_genes),'healthygene.txt')

healthy_genes<-scan(
  "healthygene.txt",
  quiet=TRUE,
  what="")

BP_healthy_genes<<-ViSEAGO::create_topGOdata(
  geneSel=healthy_genes,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_healthy_genes<-topGO::runTest(
  BP_healthy_genes,
  algorithm ="classic",
  statistic = "fisher")

BP_healthy_genes_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_healthy_genes","classic_healthy_genes")))

BP_healthy_genes_res<-as.data.frame(BP_healthy_genes_sResults@data)
write.table(BP_healthy_genes_res,"BP_healthy_genes_res.txt")


####### Boruta feature selection to shortlist key DEGs ####
set.seed(123)

vsd_df.bgc<-as.data.frame(assay(vsd_bgc))
vsd_df.bgc<-vsd_df.bgc[conditionDEG.all,]

vsd_df.bgc<-as.data.frame(t(vsd_df.bgc))
all(rownames(vsd_df.bgc)==rownames(dds@colData)) #sanity check
vsd_df.bgc$condition<-dds@colData$condition
vsd_df.bgc$condition<-as.numeric(vsd_df.bgc$condition) # change to numeric

boruta.vsd_df.bgc <- Boruta(condition~., data = vsd_df.bgc, doTrace = 2)

print(boruta.vsd_df.bgc)
#Boruta performed 99 iterations in 4.126585 mins.
#69 attributes confirmed important: ENSG00000008517, ENSG00000011028, ENSG00000011465, ENSG00000017427,
#ENSG00000026025 and 64 more;
#667 attributes confirmed unimportant: ENSG00000001626, ENSG00000007001, ENSG00000007062, ENSG00000008256,
#ENSG00000008311 and 662 more;
#84 tentative attributes left: ENSG00000013364, ENSG00000019144, ENSG00000019991, ENSG00000023445, ENSG00000025708
#and 79 more;

boruta.df1<-attStats(boruta.vsd_df.bgc)
selected.feature1<-subset(boruta.df1,boruta.df1$decision=="Confirmed")
selected.feature1$symbols<-genesymbols[rownames(selected.feature1)]
selected.feature1<-subset(selected.feature1,!selected.feature1$symbols %in% c("",'MT-TY','LINC02609','MTND1P23','MTCO3P12')) #omit some genes based on heatmap

write.csv(selected.feature1,file = "selected.feature66_finalized.csv")

final.boruta.vsd_df.bgc <- TentativeRoughFix(boruta.vsd_df.bgc) #decide for tentative features
print(final.boruta.vsd_df.bgc)
#Boruta performed 99 iterations in 4.126585 mins.
#Tentatives roughfixed over the last 99 iterations.
#122 attributes confirmed important: ENSG00000008517, ENSG00000011028, ENSG00000011465, ENSG00000013364,
#ENSG00000017427 and 117 more;
#698 attributes confirmed unimportant: ENSG00000001626, ENSG00000007001, ENSG00000007062, ENSG00000008256,
#ENSG00000008311 and 693 more;

boruta.dfFinal<-attStats(final.boruta.vsd_df.bgc)
selected.featureFinal<-subset(boruta.dfFinal,boruta.dfFinal$decision=="Confirmed")
selected.featureFinal$symbols<-genesymbols[rownames(selected.featureFinal)]
write.csv(selected.featureFinal,file = "selected.feature122_new.csv")


#heatmap of boruta genes
hm_condition.boruta<-subset(hm_mat.bgc,rownames(hm_mat.bgc) %in% rownames(selected.feature1))

hm_condition.boruta<-hm_condition.boruta[,rownames(coldat_hm)]
hm_condition.boruta<-as.data.frame(hm_condition.boruta)
hm_condition.boruta$gene<-genesymbols[rownames(hm_condition.boruta)]
hm_condition.boruta <- hm_condition.boruta %>%
  mutate(gene=case_when(!hm_condition.boruta$gene =="" ~ hm_condition.boruta$gene,
                        T ~ rownames(hm_condition.boruta)))
rownames(hm_condition.boruta)<-hm_condition.boruta$gene
hm_condition.boruta$gene<-NULL

hm_color.boruta<- viridis(100,option = "inferno")
break_hm = seq(-2, 2,length.out=100)

hm.plot.boruta_used<-pheatmap(hm_condition.boruta,scale="row",border_color = NA,color = hm_color.boruta,
                              show_rownames = T,show_colnames = F,
                              cluster_rows = T,cluster_cols =F,
                              annotation_col = coldat2_hm[,1,drop=F],
                              breaks = break_hm,
                              annotation_colors = my_color_annotation,
                              clustering_distance_rows = "correlation",
                              clustering_distance_cols = "correlation",
                              angle_col = 45) 

dev.new()
pdf("hm.plot.boruta_used.pdf",width=10,height=10)
hm.plot.boruta_used
dev.off()

selected.feature1<-selected.feature1[hm.plot.boruta_used$tree_row$order,]
write.csv(selected.feature1,file = "selected.feature66_finalized.csv")


# export NASH genes from boruta
NASHgene.boruta<-subset(hm_mat.bgc,rownames(hm_mat.bgc) %in% rownames(selected.feature1))
NASHgene.boruta<-NASHgene.boruta[,rownames(coldat2_hm)]
NASHgene.boruta<-as.data.frame(NASHgene.boruta)
NASHgene.boruta$gene<-genesymbols[rownames(NASHgene.boruta)]
NASHgene.boruta<-NASHgene.boruta[hm.plot.boruta_used$tree_row$order,]

write.csv(NASHgene.boruta[,831,drop=F],'NASHgene.boruta.csv')

#export vsd for ssGSEA
vsd.bgc_ssGSEA<-as.data.frame(assay(vsd_bgc))
vsd.bgc_ssGSEA$symbols<-genesymbols[rownames(vsd.bgc_ssGSEA)]
write.table(vsd.bgc_ssGSEA,file = "vsd.bgc_ssGSEA.txt")

#export metadata for stratification
metadata<-coldat
write.csv(metadata,file="metadata.csv")
