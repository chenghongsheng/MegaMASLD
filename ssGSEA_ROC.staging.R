setwd("./ssGSEA")

library(ssGSEA2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)

library(MASS)
library(pROC)
library(randomForest)
library(cutpointr)

###########color setting###########
color4gender<-c("blue",'pink')

pie(rep(1,4),col=brewer.pal(8,"Set3")[c(7,1,6,4)])
color4disease<-c(brewer.pal(8,"Set3")[c(7,1,6)],"#FF6347")

###########import metadata###########
metadata<-read.csv("../Deseq2_analysis/metadata.csv",header=T,check.names = F,row.names = 1)
metadata$condition<- factor(metadata$condition,levels=c("healthy",'NAFL','NASH','cirrhosis'))


###########ssGSEA using vsd###########
vsd.ssgsea<-run_ssGSEA2("./vsd.bgc.sel.gct",
                        output.prefix = "vsd.bgc.sel",
                        gene.set.databases = "./MASLDgene.gmt",
                        output.directory = "./vsd.bgc",
                        sample.norm.type = "none", 
                        weight = 0.75, 
                        correl.type = "rank", 
                        statistic = "area.under.RES",
                        output.score.type = "NES", 
                        nperm = 1000, 
                        min.overlap = 5, 
                        extended.output = TRUE, 
                        global.fdr = FALSE,
                        log.file = "./vsd.run.log")

vsd.ssgsea.score<-read.delim("./vsd.bgc/vsd.bgc.sel-scores.gct", 
                             header=F,row.names = 1,
                             comment.char="#")[-1,]

colnames(vsd.ssgsea.score) <- vsd.ssgsea.score[1,]

vsd.ssgsea.score<-as.data.frame(t(vsd.ssgsea.score))
vsd.ssgsea.score<-subset(vsd.ssgsea.score,!grepl("\\.",vsd.ssgsea.score$id))
vsd.ssgsea.score<-vsd.ssgsea.score[,-1,drop=F]
vsd.ssgsea.score$condition<-metadata$condition
vsd.ssgsea.score$`MASLD-related genes`<-as.numeric(vsd.ssgsea.score$`MASLD-related genes`)

all(rownames(vsd.ssgsea.score)==rownames(metadata)) #sanity check

vsd.ssgsea.score<-vsd.ssgsea.score[order(vsd.ssgsea.score$`MASLD-related genes`,decreasing = F),]


#heatmap
hm_mat.vsd<-t(vsd.ssgsea.score[,-2,drop=F])


hm_ssgsea.vsd<-pheatmap(hm_mat.vsd,scale="row",border_color = NA,color = hm_color,
                        show_rownames = T,show_colnames = F,
                        cluster_rows = F,cluster_cols =F,
                        annotation_col = coldat_hm.vsd,
                        #annotation_row = rowdat_hm,
                        breaks = break_hm,
                        annotation_colors = my_color_annotation,
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        angle_col = 45, cellwidth = 1.5,cellheight = 20) 

dev.new()
pdf('hm_ssgsea.vsd.pdf',width = 25)
hm_ssgsea.vsd
dev.off()

###########dotplot###########
colnames(count.ssgsea.score)<-c('score','condition')

vln_count<-ggplot(count.ssgsea.score,aes(x=condition,y=score,fill=condition))+
  geom_violin()+
  #geom_dotplot(stackdir = "center",binaxis = "y",binwidth = 0.12)+
  scale_fill_manual(values=color4disease)+
  stat_summary(fun=median,geom="crossbar",width=0.5)

dev.new()
pdf('vln_count.pdf')
vln_count
dev.off()

colnames(vsd.ssgsea.score)<-c('score','condition')

vln_vsd<-ggplot(vsd.ssgsea.score,aes(x=condition,y=score,fill=condition))+
  geom_violin()+
  #geom_dotplot(stackdir = "center",binaxis = "y",binwidth = 0.12)+
  scale_fill_manual(values=color4disease)+
  stat_summary(fun=mean,geom="crossbar",width=0.5)

dev.new()
pdf('vln_vsd.pdf')
vln_vsd
dev.off()

pairwise.t.test(vsd.ssgsea.score$score,vsd.ssgsea.score$condition,
                p.adjust.method = "holm")
#  healthy NAFL    NASH   
#  NAFL      6.6e-15 -       -      
#  NASH      < 2e-16 < 2e-16 -      
#  cirrhosis < 2e-16 < 2e-16 2.4e-09

###########ROC curve###########

data <- vsd.ssgsea.score # reads the dataset

set.seed(123)
ind <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7, 0.3))
data_validate <- data[ind, ]
data_train <- data[!ind, ]

###########ROC curve_training data###########

# Install and load necessary packages
# install.packages("MASS")
library(MASS)

# Assuming your data is in a data frame named "mydata" with columns "outcome" (ordinal) and "predictor1", "predictor2", etc. (continuous predictors)
# Adjust the column names accordingly

# Load your data
# For example:
# mydata <- read.csv("your_data.csv")

# Fit ordinal logistic regression model
model <- polr(condition ~ score, data = data_train, Hess = TRUE)

# Print summary of the model
summary(model)

# Make predictions on new data
# newdata <- data.frame(predictor1 = c(1, 2, 3), predictor2 = c(4, 5, 6))
# predict(model, newdata, type = "probs")

# If you want to check the proportional odds assumption, you can use the "brant" function from the "brant" package
# install.packages("brant")
# library(brant)
# brant(model)

# If assumptions are violated, consider alternative models or transformations
# Further diagnostics and plotting can be performed based on your specific needs

# training data
data_train$pred.probs<-predict(model,data_train,type="probs") #prob quite bad, may need to change to binomial
data_train$pred.probs<-as.data.frame(data_train$pred.probs)
data_train$pred.probs$max<-apply(data_train$pred.probs, 1, max)  

data_train <-  data_train %>%
  mutate(pred.diagnosis=case_when(data_train$pred.probs$healthy==data_train$pred.probs$max ~ "healthy",
                                  data_train$pred.probs$NAFL==data_train$pred.probs$max ~ "NAFL",
                                  data_train$pred.probs$NASH==data_train$pred.probs$max ~ "NASH",
                                  data_train$pred.probs$cirrhosis==data_train$pred.probs$max ~ "cirrhosis",
                                  T~NA))

data_train <-data_train %>%
  mutate(pred.accuracy=case_when(condition==pred.diagnosis ~ "Correct",
                                 T~ "Wrong"))

nrow(subset(data_train,data_train$pred.accuracy=="Correct"))/nrow(data_train)*100
nrow(subset(data_train,data_train$pred.accuracy=="Correct" & 
              data_train$condition=="healthy"))/nrow(subset(data_train,data_train$condition=="healthy"))*100
nrow(subset(data_train,data_train$pred.accuracy=="Correct" & 
              data_train$condition=="NAFL"))/nrow(subset(data_train,data_train$condition=="NAFL"))*100
nrow(subset(data_train,data_train$pred.accuracy=="Correct" & 
              data_train$condition=="NASH"))/nrow(subset(data_train,data_train$condition=="NASH"))*100
nrow(subset(data_train,data_train$pred.accuracy=="Correct" & 
              data_train$condition=="cirrhosis"))/nrow(subset(data_train,data_train$condition=="cirrhosis"))*100

# accuracy overall= 64.08%
# accuracy healthy= 22.22%
# accuracy NAFL = 87.59%
# accuracy NASH = 44.62%
# accuracy cirrhosis = 12.5%

data_train$pred.accuracy<-factor(data_train$pred.accuracy,levels=c("Wrong","Correct"))

summary_train.bystage<-data_train %>%
  group_by(condition,pred.accuracy) %>%
  summarise(count.stage=n())

summary_train.bystage <- as.data.frame(summary_train.bystage)%>%
  spread(pred.accuracy,count.stage)

summary_train.bystage[5,] <-c(NA,sum(summary_train.bystage[,2]),sum(summary_train.bystage[,3]))
summary_train.bystage$overall<-c(rep(NA,4),'Overall')

summary_train.bystage <- summary_train.bystage%>%
  gather(pred.accuracy,count.stage,2:3)
summary_train.bystage$pred.accuracy<-factor(summary_train.bystage$pred.accuracy,levels=c("Wrong","Correct"))


accuracy_diag.by.stage<-ggplot(summary_train.bystage,aes(x=condition,y=count.stage,fill=pred.accuracy))+
  geom_col(position = "fill",color="black")+
  scale_fill_manual(values=c("white","black"))+
  facet_grid(~overall,scales = "free",space="free")+
  ggtitle("accuracy")

dev.new()
pdf('accuracy_diag.by.stage_train.pdf')
accuracy_diag.by.stage
dev.off()


# fit binary 
data_train <- data_train %>%
  mutate(binarystage=case_when(condition %in% c("healthy","NAFL")~0,
                               T~1))

model_binary = glm(binarystage ~ score, data_train, family="binomial")     # we use the glm()-general linear model to create an instance of model

summary(model_binary)                               # summary of the model tells us the different statistical values for our independent variables after the model is created

data_train$pred.binary_outcome<-predict(model_binary,data_train,type="response")



cp<-cutpointr(data_train, pred.binary_outcome,binarystage, 
              method = maximize_metric, metric = accuracy) 
#optimal cutpoint = 0.415635; 
#AUC=0.866531; 
#accuracy=0.791837;
#sensitivity=0.703704;
#specificity=0.835366


table(Actualvalue=data_train$binarystage,Predictedvalue=data_train$pred.binary_outcome >=0.415635) 
#Predictedvalue
#Actualvalue FALSE TRUE
#         0   137   27
#         1    24   57

roc.plot.binary_train <- pROC::roc(data_train$binarystage ~ data_train$pred.binary_outcome, plot = TRUE, print.auc = TRUE)

dev.new()
pdf("ROC_binary.pdf")
plot(cp)
dev.off()

data_train <-  data_train %>%
  mutate(pred.diagnosis.binary=case_when(data_train$pred.binary_outcome>=0.415635 ~ 1,
                                         T~0))


data_train <-data_train %>%
  mutate(pred.binary.accuracy=case_when(binarystage==pred.diagnosis.binary ~ "Correct",
                                        T~ "Wrong"))

nrow(subset(data_train,data_train$pred.binary.accuracy=="Correct"))/nrow(data_train)*100
nrow(subset(data_train,data_train$pred.binary.accuracy=="Correct" & 
              data_train$binarystage=="0"))/nrow(subset(data_train,data_train$binarystage=="0"))*100
nrow(subset(data_train,data_train$pred.binary.accuracy=="Correct" & 
              data_train$binarystage=="1"))/nrow(subset(data_train,data_train$binarystage=="1"))*100


#overall accuracy = 79.18367%
#early stage/low-risk accuracy = 83.53659%
#advanced stage/low-risk accuracy = 70.37037%

data_train$pred.binary.accuracy<-factor(data_train$pred.binary.accuracy,levels=c("Wrong","Correct"))

summary_train.bybinary<-data_train %>%
  group_by(binarystage,pred.binary.accuracy) %>%
  summarise(count.binary=n())

summary_train.bybinary <- as.data.frame(summary_train.bybinary)%>%
  spread(pred.binary.accuracy,count.binary)

summary_train.bybinary[3,] <-c(NA,sum(summary_train.bybinary[,2]),sum(summary_train.bybinary[,3]))
summary_train.bybinary$overall<-c(rep(NA,2),'Overall')

summary_train.bybinary <- summary_train.bybinary%>%
  gather(pred.binary.accuracy,count.binary,2:3)
summary_train.bybinary$pred.binary.accuracy<-factor(summary_train.bybinary$pred.binary.accuracy,levels=c("Wrong","Correct"))
summary_train.bybinary$binarystage<-as.factor(summary_train.bybinary$binarystage)

accuracy_diag.by.binary<-ggplot(summary_train.bybinary,aes(x=binarystage,y=count.binary,fill=pred.binary.accuracy))+
  geom_col(position = "fill",color="black")+
  scale_fill_manual(values=c("white","black"))+
  facet_grid(~overall,scales = "free",space="free")+
  ggtitle("accuracy")



dev.new()
pdf('accuracy_diag.by.binarystage_train.pdf')
accuracy_diag.by.binary
dev.off()


###########ROC curve_validation data###########
data_validate$pred.probs<-predict(model,data_validate,type="probs")
data_validate$pred.probs<-as.data.frame(data_validate$pred.probs)
data_validate$pred.probs$max<-apply(data_validate$pred.probs, 1, max)  

data_validate <-  data_validate %>%
  mutate(pred.diagnosis=case_when(data_validate$pred.probs$healthy==data_validate$pred.probs$max ~ "healthy",
                                  data_validate$pred.probs$NAFL==data_validate$pred.probs$max ~ "NAFL",
                                  data_validate$pred.probs$NASH==data_validate$pred.probs$max ~ "NASH",
                                  data_validate$pred.probs$cirrhosis==data_validate$pred.probs$max ~ "cirrhosis",
                                  T~NA))

data_validate <-data_validate %>%
  mutate(pred.accuracy=case_when(condition==pred.diagnosis ~ "Correct",
                                 T~ "Wrong"))

nrow(subset(data_validate,data_validate$pred.accuracy=="Correct"))/nrow(data_validate)*100
nrow(subset(data_validate,data_validate$pred.accuracy=="Correct" & 
              data_validate$condition=="healthy"))/nrow(subset(data_validate,data_validate$condition=="healthy"))*100
nrow(subset(data_validate,data_validate$pred.accuracy=="Correct" & 
              data_validate$condition=="NAFL"))/nrow(subset(data_validate,data_validate$condition=="NAFL"))*100
nrow(subset(data_validate,data_validate$pred.accuracy=="Correct" & 
              data_validate$condition=="NASH"))/nrow(subset(data_validate,data_validate$condition=="NASH"))*100
nrow(subset(data_validate,data_validate$pred.accuracy=="Correct" & 
              data_validate$condition=="cirrhosis"))/nrow(subset(data_validate,data_validate$condition=="cirrhosis"))*100

# accuracy overall= 57.4359%
# accuracy healthy= 6.557377%
# accuracy NAFL = 82.7044%
# accuracy NASH = 35.67251%
# accuracy cirrhosis = 22.85714%

data_validate$pred.accuracy<-factor(data_validate$pred.accuracy,levels=c("Wrong","Correct"))

summary_validate.bystage<-data_validate %>%
  group_by(condition,pred.accuracy) %>%
  summarise(count.stage=n())

summary_validate.bystage <- as.data.frame(summary_validate.bystage)%>%
  spread(pred.accuracy,count.stage)

summary_validate.bystage[5,] <-c(NA,sum(summary_validate.bystage[,2]),sum(summary_validate.bystage[,3]))
summary_validate.bystage$overall<-c(rep(NA,4),'Overall')

summary_validate.bystage <- summary_validate.bystage%>%
  gather(pred.accuracy,count.stage,2:3)
summary_validate.bystage$pred.accuracy<-factor(summary_validate.bystage$pred.accuracy,levels=c("Wrong","Correct"))


accuracy_diag.by.stage<-ggplot(summary_validate.bystage,aes(x=condition,y=count.stage,fill=pred.accuracy))+
  geom_col(position = "fill",color="black")+
  scale_fill_manual(values=c("white","black"))+
  facet_grid(~overall,scales = "free",space="free")+
  ggtitle("accuracy")

dev.new()
pdf('accuracy_diag.by.stage_validate.pdf')
accuracy_diag.by.stage
dev.off()


#binary
data_validate <- data_validate %>%
  mutate(binarystage=case_when(condition %in% c("healthy","NAFL")~0,
                               T~1))

data_validate$pred.binary_outcome<-predict(model_binary,data_validate,type="response")

cp_validate<-cutpointr(data_validate, pred.binary_outcome,binarystage, 
                       method = maximize_metric, metric = accuracy)

#optimal cutpoint = 0.415635; follow the train model 
#AUC=0.7815; 
#accuracy=;0.7470085
#sensitivity=0.5970874;
#specificity=0.8284960

dev.new()
pdf("ROC_binary_validate.pdf")
plot(cp_validate)
dev.off()


table(Actualvalue=data_validate$binarystage,Predictedvalue=data_validate$pred.binary_outcome >=0.415635) 
#Predictedvalue
#Actualvalue FALSE TRUE
#         0   314   65
#         1    83   123

roc.plot.binary_validate <- pROC::roc(data_validate$binarystage ~ data_validate$pred.binary_outcome, plot = TRUE, print.auc = TRUE)

dev.new()
pdf("ROC_binary.pdf")
plot(cp_validate)
dev.off()

data_validate <-  data_validate %>%
  mutate(pred.diagnosis.binary=case_when(data_validate$pred.binary_outcome>=0.415635 ~ 1,
                                         T~0))


data_validate <-data_validate %>%
  mutate(pred.binary.accuracy=case_when(binarystage==pred.diagnosis.binary ~ "Correct",
                                        T~ "Wrong"))

nrow(subset(data_validate,data_validate$pred.binary.accuracy=="Correct"))/nrow(data_validate)*100
nrow(subset(data_validate,data_validate$pred.binary.accuracy=="Correct" & 
              data_validate$binarystage=="0"))/nrow(subset(data_validate,data_validate$binarystage=="0"))*100
nrow(subset(data_validate,data_validate$pred.binary.accuracy=="Correct" & 
              data_validate$binarystage=="1"))/nrow(subset(data_validate,data_validate$binarystage=="1"))*100
#overall accuracy = 74.70085%
#early stage/low-risk accuracy = 82.8496%
#advanced stage/low-risk accuracy = 59.70874%

data_validate$pred.binary.accuracy<-factor(data_validate$pred.binary.accuracy,levels=c("Wrong","Correct"))

summary_validate.bybinary<-data_validate %>%
  group_by(binarystage,pred.binary.accuracy) %>%
  summarise(count.binary=n())

summary_validate.bybinary <- as.data.frame(summary_validate.bybinary)%>%
  spread(pred.binary.accuracy,count.binary)

summary_validate.bybinary[3,] <-c(NA,sum(summary_validate.bybinary[,2]),sum(summary_validate.bybinary[,3]))
summary_validate.bybinary$overall<-c(rep(NA,2),'Overall')

summary_validate.bybinary <- summary_validate.bybinary%>%
  gather(pred.binary.accuracy,count.binary,2:3)
summary_validate.bybinary$pred.binary.accuracy<-factor(summary_validate.bybinary$pred.binary.accuracy,levels=c("Wrong","Correct"))
summary_validate.bybinary$binarystage<-as.factor(summary_validate.bybinary$binarystage)

accuracy_diag.by.binary<-ggplot(summary_validate.bybinary,aes(x=binarystage,y=count.binary,fill=pred.binary.accuracy))+
  geom_col(position = "fill",color="black")+
  scale_fill_manual(values=c("white","black"))+
  facet_grid(~overall,scales = "free",space="free")+
  ggtitle("accuracy")

dev.new()
pdf('accuracy_diag.by.binarystage_validate.pdf')
accuracy_diag.by.binary
dev.off()
