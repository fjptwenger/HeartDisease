rm(list=ls())

#----------PART1 GTEx data---------
#read phenotype information,and extract COAD data
#download website: https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz
GTEx_pheno<-read.table('GTEX_phenotype',sep='\t',header = T)
GTEx_pheno_colon<-GTEx_pheno[substr(GTEx_pheno$body_site_detail..SMTSD.,1,5)=='Colon',]
rownames(GTEx_pheno_colon)<-GTEx_pheno_colon$Sample

options(stringsAsFactors = F)
library(data.table)
#download website:https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_gene_expected_count.gz
GTEx_data<-fread('gtex_gene_expected_count',sep='\t',header = T)
GTEx_data_new<-as.data.frame(GTEx_data)
GTEx_data_colon<-GTEx_data_new[,colnames(GTEx_data_new)  %in% rownames(GTEx_pheno_colon)]
rownames(GTEx_data_colon)<-GTEx_data_new$sample
GTEx_data_colon$id<-GTEx_data_new$sample

GTEx_pheno_colon_select<-GTEx_pheno_colon[GTEx_pheno_colon$Sample %in% colnames(GTEx_data_new),]
#ID transfer, download website:https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap
expr_probe<-read.table('gencode.v23.annotation.gene.probemap',sep='\t',header = T)
probe2symbol<-expr_probe[,1:2]

GTEx_data_colon_gene<-merge(GTEx_data_colon,probe2symbol,by='id')
GTEx_data_colon_gene$id<-NULL
GTEx_data_colon_gene_dedup <- aggregate(x = GTEx_data_colon_gene,by = list(GTEx_data_colon_gene$gene), FUN = max)
GTEx_data_colon_gene_dedup$Group.1<-NULL
rownames(GTEx_data_colon_gene_dedup)<-GTEx_data_colon_gene_dedup$gene

library(dplyr)
GTEx_data_colon_gene_dedup<-GTEx_data_colon_gene_dedup %>% na.omit()

#If the read counts of a gene in more than 50% of the samples is less than 10, it will be filtered out.
qualified_genes<-c()
for (gene_in_sheet in rownames(GTEx_data_colon_gene_dedup)){
  qualification<-GTEx_data_colon_gene_dedup[gene_in_sheet,]<(log(10,2)+1)
  if (sum(qualification)<0.5*length(GTEx_data_colon_gene_dedup)){
    qualified_genes<-append(qualified_genes,gene_in_sheet)
  }
}
GTEx_data_colon_filt<-GTEx_data_colon_gene_dedup[qualified_genes,]

#----------PART2 TCGA data---------
#1. expression data of COAD
#phenotype information
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.GDC_phenotype.tsv.gz
TCGA_colon_pheno<-read.table('TCGA-COAD.GDC_phenotype.tsv',header = T,sep='\t',quote = '"')

#RNA-seq read
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.htseq_counts.tsv.gz
TCGA_colon_expr<-read.table('TCGA-COAD.htseq_counts.tsv',header=T,sep='\t')
#filt out type B or type C samples
TCGA_colon_expr_sample<-TCGA_colon_expr[,substr(colnames(TCGA_colon_expr),16,16)!='B' & substr(colnames(TCGA_colon_expr),16,16)!='C']
TCGA_colon_expr_sample$id<-TCGA_colon_expr$Ensembl_ID
TCGA_colon_expr_sample$Ensembl_ID<-NULL

#ID transfer 
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap
expr_probe<-read.table('gencode.v22.annotation.gene.probeMap',sep='\t',header = T)
probe2symbol<-expr_probe[,1:2]
TCGA_colon_expr_gene<-merge(TCGA_colon_expr_sample,probe2symbol,by='id')
TCGA_colon_expr_gene$id<-NULL
TCGA_colon_expr_gene_dedup <- aggregate(x = TCGA_colon_expr_gene,by = list(TCGA_colon_expr_gene$gene), FUN = max)
TCGA_colon_expr_gene_dedup$Group.1<-NULL
rownames(TCGA_colon_expr_gene_dedup)<-TCGA_colon_expr_gene_dedup$gene

colnames(TCGA_colon_expr_gene_dedup)<-gsub("\\.","-",colnames(TCGA_colon_expr_gene_dedup))
colnames(TCGA_colon_expr_gene_dedup)<-substr(colnames(TCGA_colon_expr_gene_dedup),1,15)
TCGA_colon_expr_gene_dedup$gene<-NULL

#If the read counts of a gene in more than 50% of the samples is less than 10, it will be filtered out.
qualified_genes<-c()
for (gene_in_sheet in rownames(TCGA_colon_expr_gene_dedup)){
  qualification<-TCGA_colon_expr_gene_dedup[gene_in_sheet,]<(log(10,2)+1)
  if (sum(qualification)<0.5*length(TCGA_colon_expr_gene_dedup)){
    qualified_genes<-append(qualified_genes,gene_in_sheet)
  }
}
TCGA_colon_expr_filt<-TCGA_colon_expr_gene_dedup[qualified_genes,]

group_list<-ifelse(substr(colnames(TCGA_colon_expr_filt),14,15)<10,'Tumor','Normal')

#2. expression data--READ
#read phenotype information
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.GDC_phenotype.tsv.gz
TCGA_READ_pheno<-read.table('TCGA-READ.GDC_phenotype.tsv',header = T,sep='\t',quote = '"')

#read RNA-seq reads
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.htseq_counts.tsv.gz
TCGA_READ_expr<-read.table('TCGA-READ.htseq_counts.tsv',header=T,sep='\t')
TCGA_READ_expr$id<-TCGA_READ_expr$Ensembl_ID
#ID transfer
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap
expr_probe<-read.table('rawdata/gencode.v22.annotation.gene.probeMap',sep='\t',header = T)
probe2symbol<-expr_probe[,1:2]
TCGA_READ_expr_gene<-merge(TCGA_READ_expr,probe2symbol,by='id')
TCGA_READ_expr_gene$id<-NULL
TCGA_READ_expr_gene$Ensembl_ID<-NULL
TCGA_READ_expr_gene_dedup <- aggregate(x = TCGA_READ_expr_gene,by = list(TCGA_READ_expr_gene$gene), FUN = max)
TCGA_READ_expr_gene_dedup$Group.1<-NULL
rownames(TCGA_READ_expr_gene_dedup)<-TCGA_READ_expr_gene_dedup$gene
#If the read counts of a gene in more than 50% of the samples is less than 10, it will be filtered out.
qualified_genes<-c()
for (gene_in_sheet in rownames(TCGA_READ_expr_gene_dedup)){
  qualification<-TCGA_READ_expr_gene_dedup[gene_in_sheet,]<(log(10,2)+1)
  if (sum(qualification)<0.5*length(TCGA_READ_expr_gene_dedup)){
    qualified_genes<-append(qualified_genes,gene_in_sheet)
  }
}
TCGA_READ_expr_filt<-TCGA_READ_expr_gene_dedup[qualified_genes,]
colnames(TCGA_READ_expr_filt)<-gsub("\\.","-",colnames(TCGA_READ_expr_filt))
TCGA_READ_expr_filt$gene<-NULL
colnames(TCGA_READ_expr_filt)<-substr(colnames(TCGA_READ_expr_filt),1,15)
group_list<-ifelse(substr(colnames(TCGA_READ_expr_filt),14,15)<10,'Tumor','Normal')

#merge phenotype information
TCGA_pheno_data<-rbind(TCGA_colon_pheno[,c('submitter_id.samples','age_at_initial_pathologic_diagnosis','gender.demographic','tumor_stage.diagnoses','pathologic_T','pathologic_M','pathologic_N')],  
                      TCGA_READ_pheno[,c('submitter_id.samples','age_at_initial_pathologic_diagnosis','gender.demographic','tumor_stage.diagnoses','pathologic_T','pathologic_M','pathologic_N')]
                      )

TCGA_pheno_data_filt<-TCGA_pheno_data[substr(TCGA_pheno_data$submitter_id.samples,16,16)=='A',]
rownames(TCGA_pheno_data_filt)<-substr(TCGA_pheno_data_filt$submitter_id.samples,1,15)
TCGA_pheno_data_filt2<-TCGA_pheno_data_filt[colnames(COAD_READ_expr),]
colnames(TCGA_pheno_data_filt2)=c('sampleid','Age','Gender','Stage','T','M','N')
library(stringr)
TCGA_pheno_data_filt2 <- TCGA_pheno_data_filt2 %>% na.omit()
TCGA_pheno_data_filt2$Age<-factor(ifelse(TCGA_pheno_data_filt2$Age >60,'>60','<=60'))
TCGA_pheno_data_filt2$Group<-ifelse(substr(rownames(TCGA_pheno_data_filt2),14,15)<10,'Tumor','Normal')
TCGA_pheno_data_filt2$Stage<-factor(toupper(str_extract(TCGA_pheno_data_filt2$Stage, "i+v*")))
TCGA_pheno_data_filt2$T<-factor(str_extract(TCGA_pheno_data_filt2$T, "T[0-9]+"))
TCGA_pheno_data_filt2$M <- factor(str_extract(TCGA_pheno_data_filt2$M, "M[0-9]*X*"))
TCGA_pheno_data_filt2$N<-factor(str_extract(TCGA_pheno_data_filt2$N, "N[0-9]"))

TCGA_pheno_data_filt2$Group<-factor(TCGA_pheno_data_filt2$Group,
                levels=c('Tumor','Normal'),
                labels=c("Tumor", # 第一个作为参考组
                         "Normal"))

library(table1)
table1(~ factor(Gender) + factor(Age) + factor(Stage) +
         factor(T) +  factor(M) +  factor(N) | Group, data=TCGA_pheno_data_filt2)

#merge COAD and READ reads
commongene<-intersect(rownames(TCGA_colon_expr_filt),rownames(TCGA_READ_expr_filt))
COAD_READ_expr<-cbind(TCGA_colon_expr_filt[commongene,],TCGA_READ_expr_filt[commongene,])

#PCA of merged COAD&READ expression data ===Figure S1===
group<-ifelse(colnames(COAD_READ_expr) %in% colnames(TCGA_colon_expr_filt),ifelse(substr(colnames(COAD_READ_expr),14,15)<10,'COAD-Tumor','COAD-Normal'),ifelse(substr(colnames(COAD_READ_expr),14,15)<10,'READ-Tumor','READ-Normal'))
df=as.data.frame(t(COAD_READ_expr))#行列转换	  
df$group=group
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = "group")+theme_bw()#画图


#DEG analysis
source('exp_diff_limma.R')
DEG<-processData(COAD_READ_expr)

#画火山图
ps<-DEG$P.Value
rs<-DEG$logFC
cs<-rep('darkgray',length(ps))
cs[ps<0.05 & rs>0.585]<-'red'
cs[ps<0.05 & rs<(-0.585)]<-'green'

plot(rs,-log10(ps),col=cs,pch=19,cex=0.5,xlab = 'log2FC',ylab = '-log10(P-value)')


#3. methylation data
#COAD raw data
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.methylation450.tsv.gz
TCGA_COAD_methy<-read.table('COAD_Methylation450.xls',header = T,sep='\t')

#READ raw data
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.methylation450.tsv.gz
TCGA_READ_methy<-read.table('rawdata/TCGA-READ.methylation450.tsv',header = T,sep='\t')
TCGA_READ_methy$sample<-TCGA_READ_methy$Composite.Element.REF
TCGA_READ_methy$Composite.Element.REF<-NULL

#merge methylation data
TCGA_methy_expr<-merge(TCGA_COAD_methy,TCGA_READ_methy,by='sample')
colnames(TCGA_methy_expr)<-gsub('\\.','-',colnames(TCGA_methy_expr))
colnames(TCGA_methy_expr)<-substr(colnames(TCGA_methy_expr),1,15)
suppressMessages(library(ChAMP))
suppressMessages(library(stringr))
rownames(TCGA_methy_expr)=TCGA_methy_expr[,1]
TCGA_methy_expr<-TCGA_methy_expr[,-1]
library(dplyr) 
TCGA_methy_expr_filt<- TCGA_methy_expr  %>% na.omit()

group_list<-ifelse(as.numeric(substr(colnames(TCGA_methy_expr_filt),14,15))<10,"Tumor","Normal")
TCGA_methy_pd<-data.frame(colnames(TCGA_methy_expr_filt),substr(colnames(TCGA_methy_expr_filt),1,12))
colnames(TCGA_methy_pd)<-c("sampleID",'subdivision')
TCGA_methy_pd$group_list<-group_list
TCGA_beta=as.matrix(TCGA_methy_expr_filt)
memory.limit(2000000) 
TCGA_beta_filter<-champ.filter(beta = TCGA_beta ,pd = TCGA_methy_pd) #This step has been automatically filtered raw data

#Methylation annotation and difference analysis
source("methy_diff.R")
TCGA_DMP_siggene<-processMethy(beta=TCGA_beta_filter,logfc_cutoff=0.3)

#draw the phenotype information table
TCGA_methy_pheno<-TCGA_pheno_data_filt[TCGA_methy_pd$sampleID,] 
colnames(TCGA_methy_pheno)=c('sampleid','Age','Gender','Stage','T','M','N')
library(stringr)
TCGA_methy_pheno <- TCGA_methy_pheno %>% na.omit()
TCGA_methy_pheno$Age<-factor(ifelse(TCGA_methy_pheno$Age >60,'>60','<=60'))
TCGA_methy_pheno$Group<-ifelse(substr(rownames(TCGA_methy_pheno),14,15)<10,'Tumor','Normal')
TCGA_methy_pheno$Stage<-factor(toupper(str_extract(TCGA_methy_pheno$Stage, "i+v*")))
TCGA_methy_pheno$T<-factor(str_extract(TCGA_methy_pheno$T, "T[0-9]+"))
TCGA_methy_pheno$M <- factor(str_extract(TCGA_methy_pheno$M, "M[0-9]*X*"))
TCGA_methy_pheno$N<-factor(str_extract(TCGA_methy_pheno$N, "N[0-9]"))

TCGA_methy_pheno$Group<-factor(TCGA_methy_pheno$Group,
                               levels=c('Tumor','Normal'),
                               labels=c("Tumor", # 第一个作为参考组
                                        "Normal"))

library(table1)
table1(~ factor(Gender) + factor(Age) + factor(Stage) +
         factor(T) +  factor(M) +  factor(N) | Group, data=TCGA_methy_pheno)

#PART 3. methylation driven genes analysis
dgeSigAll<-DEG[DEG$change!='NOT',]
common_gene<-intersect(rownames(dgeSigAll),TCGA_DMP_siggene$gene)

TCGA_common_sample<-intersect(colnames(COAD_READ_expr),colnames(TCGA_methy_expr_filt))
TCGA_tumor_sample<-TCGA_common_sample[substr(TCGA_common_sample,14,15)<10]
TCGA_normal_sample<-setdiff(TCGA_common_sample,TCGA_tumor_sample)

GE_cancer<-COAD_READ_expr[common_gene,TCGA_tumor_sample]

cg2gene<-data.frame(rownames(TCGA_DMP_siggene),TCGA_DMP_siggene$gene)
colnames(cg2gene)=c('cg','GeneName')
methy_data<-TCGA_beta_filter$beta
methy_gene_data<-methy_data[rownames(methy_data) %in% rownames(TCGA_DMP_siggene),]
methy_gene_add<-data.frame(rownames(methy_gene_data),methy_gene_data)
methy_gene_add$cg<-methy_gene_add[,1]
methy_gene_add<-methy_gene_add[,-1]
methy_gene_merge<-merge(cg2gene, methy_gene_add, by = "cg")
methy_gene_merge$cg<-NULL
methy_Gene_driver<-methy_gene_merge[methy_gene_merge$GeneName %in% common_gene,]
methy_Gene_driver_dedup <- aggregate(x = as.matrix(methy_Gene_driver),by = list(methy_Gene_driver$GeneName), FUN = max)
rownames(methy_Gene_driver_dedup)<-methy_Gene_driver_dedup[,1]
methy_Gene_driver_dedup$Group.1<-NULL
colnames(methy_Gene_driver_dedup)<-gsub('\\.','-',colnames(methy_Gene_driver_dedup))
methy_Gene_driver_dedup$GeneName<-NULL
MET_cancer<-methy_Gene_driver_dedup[common_gene,TCGA_tumor_sample]
MET_normal<-methy_Gene_driver_dedup[common_gene,TCGA_tumor_sample]

MET_cancer_new<-apply(MET_cancer, 2, as.numeric)
rownames(MET_cancer_new)<-rownames(MET_cancer)
MET_cancer_new<-as.matrix(MET_cancer_new)

GE_cancer_new<-apply(GE_cancer, 2, as.numeric)
rownames(GE_cancer_new)<-rownames(GE_cancer)
GE_cancer_new<-as.matrix(GE_cancer_new)

MET_normal_new<-apply(MET_normal, 2, as.numeric)
rownames(MET_normal_new)<-rownames(MET_normal)
MET_normal_new<-as.matrix(MET_normal_new)

source("doMethylMix.R")
MET_result<-processMethylMix(MET_cancer_new, GE_cancer_new, MET_normal_new)

driver_genes<-MET_result$MethylationDrivers
#extract driver genes information
TCGA_DMP_siggene_driver<-TCGA_DMP_siggene[TCGA_DMP_siggene$gene %in% driver_genes,c('gene','logFC','P.Value','adj.P.Val','Tumor_AVG','Normal_AVG')]
TCGA_DMP_siggene_driver$gene<-as.character(TCGA_DMP_siggene_driver$gene)
TCGA_DMP_siggene_driver_dedup<-aggregate(x = TCGA_DMP_siggene_driver,by = list(TCGA_DMP_siggene_driver$gene), FUN = max)   
TCGA_DMP_siggene_driver_dedup$Group.1<-NULL

#Combining GTEX data, using dirver gene to build prediction model
GTEx_data_colon_drivergene<-GTEx_data_colon_gene_dedup[driver_genes,] %>% na.omit
GTEx_data_colon_drivergene$gene<-NULL
COAD_READ_expr_drivergene<-COAD_READ_expr[driver_genes,]
COAD_READ_expr_drivergene$gene<-NULL
merge_exprs<-cbind(GTEx_data_colon_drivergene,COAD_READ_expr_drivergene)

#PCA of merged TCGA& GTEx expression data ===Figure S1===
group<-c(rep('GTEx_COAD_normal',length(colnames(GTEx_data_colon_drivergene))),ifelse(as.numeric(substr(colnames(COAD_READ_expr_drivergene),14,15))<10,'TCGA_Tumor','TCGA_Normal'))
df=as.data.frame(t(merge_exprs))	  
df$group=group
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = "group")+theme_bw()


#draw heatmap using drivergene
suppressMessages(library(pheatmap))
plot_matrix<-merge_exprs
group_list<-c(rep('GTEx_Normal',length(colnames(GTEx_data_colon_drivergene))),ifelse(as.numeric(substr(colnames(COAD_READ_expr_drivergene),14,15))<10,'TCGA_Tumor','TCGA_Normal'))
plot_matrix<-scale(plot_matrix)#画热图要先归一化 
annotation_col <- data.frame(Group=group_list) 
rownames(annotation_col) <- colnames(plot_matrix)
ann_colors = list(Group = c(TCGA_Tumor="#eb4025", TCGA_Normal="#bb51f7",GTEx_Normal="#78f47b"))
tiff(file='drivergene_heatmap.tiff',width=15,height=7,units="cm",compression="lzw",bg="white",res=500)
p=pheatmap(plot_matrix,show_colnames = F,show_rownames=T,
         annotation_col = annotation_col,
         border_color=NA,
         color =  colorRampPalette(c("green","black", "red"), bias = 1.2)(100),
         annotation_colors = ann_colors)
print(p)
dev.off()


#Correlation between driver genes
cor_matrix<-t(merge_exprs)
factor_Corr <- cor(cor_matrix)
library(corrplot)
tiff(file='drivergene_cor.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=500)
corrplot(factor_Corr,method="number")
dev.off()

targetgeneinfo<-dgeSigAll[driver_genes,]
group_list<-c(rep('Normal',length(colnames(GTEx_data_colon_drivergene))),ifelse(as.numeric(substr(colnames(COAD_READ_expr_drivergene),14,15))<10,'Tumor','Normal'))
driver_gene_matrix<-data.frame(group_list,t(merge_exprs))

gene_in=array()
for (gene in driver_genes){
  uni_gene_in<-data.frame(gene,driver_gene_matrix[,gene],rownames(driver_gene_matrix),group_list)
  gene_in<-rbind(gene_in,uni_gene_in)
}

colnames(gene_in)<-c('gene','expression','sample','group')

gene_in<-gene_in %>% na.omit()

#draw MDGs 
library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
p <- ggboxplot(gene_in, x = "gene", y = "expression",
               color = 'group', palette = "nejm",
               add = "jitter")
p + stat_compare_means(aes(group = group),method = "t.test",label = "p.signif")+rotate_x_text(angle = 45)



learn_matrix<-driver_gene_matrix
##T-test
pvalue<-sapply(2:ncol(learn_matrix),function(i){t.test(learn_matrix[which(group_list=='Tumor'),i],learn_matrix[which(group_list=='Normal'),i])[[3]]})
#The p value of t test was corrected to reduce false positive
qvlaue<-p.adjust(pvalue)
##select gene with P<0.05
learn_matrix_filt<-learn_matrix[,which(qvlaue<0.05)+1]
learn_matrix_filt<-as.data.frame(learn_matrix_filt)
#order by group_list
learn_matrix_filt<-learn_matrix_filt[rownames(learn_matrix),]
learn_matrix_new<-cbind(group_list,learn_matrix_filt)

data<-learn_matrix_new
data$group_list<-as.factor(group_list)
save(data,file = "dataforprediction.Rdata")

load('dataforprediction.Rdata')
#build prediction model using RF algorithm
#BiocManager::install('randomForest')
#BiocManager::install('caret')
library(randomForest)
library(caret)
library(pROC)
inTrain<-createDataPartition(y=data[,1],p=0.25,list = F)
test<-data[inTrain,] 
train<-data[-inTrain,]

#MDGs: "ARHGAP20" "CLDN1"    "KCNJ12"   "STK33"    "EPHX4"    "SPTBN5"   "KRT7"     "FAM150A"  "TCN1"     "LY6G6D"
#calculate ROC of each MDG
for(i in driver_genes){
  formula<-as.formula(paste0('group_list~',i,sep=''))
  model<-randomForest(formula,data=train,proximity=T,importance=T)
  model_pre<-predict(model,newdata = test,type="prob")
  rf_roc <- roc(test$group_list,model_pre[,1])
  tiff(file=paste('ROC',paste0(i,'.tiff',sep=''),sep='/'),width=12,height=12,units="cm",compression="lzw",bg="white",res=600)
  p=plot(rf_roc, print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main=i)
  print(p)
  dev.off()
}

#calculate ROC of MDGs panel
formula<-as.formula(paste0('group_list~',paste(colnames(learn_matrix_filt),sep='',collapse = '+')))
model<-randomForest(formula,data=train,proximity=T,importance=T)
model_pre<-predict(model,newdata = test,type="prob")
rf_roc <- roc(test$group_list,model_pre[,1])
tiff(file='overall_ROC.tiff',width=12,height=12,units="cm",compression="lzw",bg="white",res=600)
p=plot(rf_roc, print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='10 gene panel')
print(p)
dev.off()

#drawing heatmap 
#traing data
train_predict<-as.matrix(model$predicted)
colnames(train_predict)<-c('predict_status')
train_real<-train[rownames(train_predict),]
train_predict_real<-cbind(train_predict,train_real)
train_predict_real$real_status<-train_predict_real$group_list
train_predict_real$group_list<-NULL

#testing data
test_predict<-data.frame(rownames(model_pre),ifelse(model_pre[,1]<0.5,'Pre_tumor','Pre_normal'))
colnames(test_predict)<-c('samples','predict_status')
test_real<-test[rownames(test_predict),]
test_predict_real<-cbind(test_predict,test_real)
test_predict_real$real_status<-test_predict_real$group_list
test_predict_real$group_list<-NULL
test_predict_real$samples<-NULL

library(pheatmap)
#plot_matrix<-train_predict_real
#plot_matrix$predict_status<-ifelse(train_predict_real$predict_status=='Normal','Pre_normal','Pre_cancer')
plot_matrix<-test_predict_real
plot_matrix$predict_status<-ifelse(test_predict_real$predict_status=='Normal','Pre_normal','Pre_cancer')

predict_status<-c('navy','skyblue')
real_status<-c('red','green')
ann_colors=list(predict_status=predict_status,real_status=real_status)
pcol.map<-function(s){if(s=='Pre_normal') 1 else 2}
rcol.map<-function(s){if(s=='Normal') 3 else 4}
pcolors<-unlist(lapply(plot_matrix$predict_status, pcol.map))
rcolors<-unlist(lapply(plot_matrix$real_status, rcol.map))

annotation<-data.frame(plot_matrix$real_status,plot_matrix$predict_status)
rownames(annotation)<-rownames(plot_matrix)
colnames(annotation)<-c('Real status','Predict status')
tiff(file='heatmaptest.tiff',width=15,height=6,units="cm",compression="lzw",bg="white",res=600)
p=pheatmap(scale(t(plot_matrix[,2:11])),show_colnames = F,show_rownames=T,
         annotation = annotation,
         fontsize_row = 6,
         fontsize = 6,
         border_color=NA,
         color =  colorRampPalette(c("green","black", "red"), bias = 1.2)(100),
         annotation_colors = ann_colors)
print(p)
dev.off()

#calculate the importance of each gene
varImpPlot(model, main = "variable importance",pch=15)

#GSE39582 for validation the model
library(GEOquery)
f='GSE39582_eSet.Rdata'

if(!file.exists(f)){
  gset <- getGEO('GSE39582', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}

table(rownames(exprSet) %in% ids$probe_id) 
exprSet = exprSet[rownames(exprSet) %in% ids$probe_id,]
tmp = by(exprSet, ids$symbol, function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
exprSet = exprSet[rownames(exprSet) %in% probes,] 
#id transfer 
ids = ids[match(rownames(exprSet),ids$probe_id),]
rownames(exprSet) <- ids$symbol
GSE39582_drivergenes<-exprSet[rownames(exprSet) %in% driver_genes,]
common_drivergene<-rownames(GSE39582_drivergenes)
GSE39582_formodel<-as.data.frame(t(GSE39582_drivergenes))
GSE39582_formodel$group_list<-ifelse(substr(clinval$source_name_ch1,18,20)=='non','Normal','Tumor')

GEO_data<-data[,common_drivergene]
GEO_data$group_list<-data$group_list

#construction of predict model
svmforGEO<-randomForest(group_list~.,data=GEO_data,proximity=T,importance=T)
model_pre<-predict(svmforGEO,newdata = GSE39582_formodel,type="prob")
rf_roc <- roc(GSE39582_formodel$group_list,model_pre[,1])
tiff(file='GSE39582_ROC.tiff',width=12,height=12,units="cm",compression="lzw",bg="white",res=600)
p=plot(rf_roc, print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='GSE39582')
print(p)
dev.off()


#draw table for GSE39582 
GEO_pheno_data<-data.frame(clinval$geo_accession,clinval$`age.at.diagnosis (year):ch1`,clinval$`Sex:ch1`,clinval$`tnm.stage:ch1`,clinval$`tnm.t:ch1`,clinval$`tnm.m:ch1`,clinval$`tnm.n:ch1`,clinval$source_name_ch1)
rownames(GEO_pheno_data)<-GEO_pheno_data$clinval.geo_accession
GEO_pheno_data<-GEO_pheno_data[rownames(GEO_pheno_data) %in% colnames(exprSet),]

colnames(GEO_pheno_data)=c('sampleid','Age','Gender','Stage','T','M','N','Group')
library(stringr)
GEO_pheno_data$Age<-factor(ifelse(GEO_pheno_data$Age >60,'>60','<=60'))
GEO_pheno_data$Group<-ifelse(substr(GEO_pheno_data$Group,18,20)=='non','Normal','Tumor')
GEO_pheno_data$Stage<-factor(toupper(str_extract(GEO_pheno_data$Stage, "[1-9]+")))
GEO_pheno_data$T<-factor(str_extract(GEO_pheno_data$T, "T[1-9]+"))
GEO_pheno_data$M <- factor(str_extract(GEO_pheno_data$M, "M[0-9]*X*"))
GEO_pheno_data$N<-factor(str_extract(GEO_pheno_data$N, "N[0-2]"))

GEO_pheno_data$Group<-factor(GEO_pheno_data$Group,
                                    levels=c('Tumor','Normal'),
                                    labels=c("Tumor", # 第一个作为参考组
                                             "Normal"))

library(table1)
table1(~ factor(Gender) + factor(Age) + factor(Stage) +
         factor(T) +  factor(M) +  factor(N) | Group, data=GEO_pheno_data)

#read survival data of COAD
#website: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.survival.tsv
COAD_survival<-read.table('TCGA-COAD.survival.tsv',header = T,sep = '\t')
COAD_survival<-COAD_survival[substr(COAD_survival$sample,16,16)=='A',]
COAD_survival$sample<-substr(COAD_survival$sample,1,15)

#read survival data of READ
#website:https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.survival.tsv
READ_survival<-read.table('TCGA-READ.survival.tsv',header = T,sep = '\t')
READ_survival<-READ_survival[substr(READ_survival$sample,16,16)=='A',]
READ_survival$sample<-substr(READ_survival$sample,1,15)

TCGA_survival<-rbind(COAD_survival,READ_survival)
TCGA_survival$X_PATIENT<-NULL
rownames(TCGA_survival)<-TCGA_survival$sample

#select the samples have survival information
common_sample<-intersect(rownames(TCGA_expr),rownames(TCGA_survival))
survival_gene_data<-cbind(TCGA_expr[common_sample,],TCGA_survival[common_sample,])
survival_gene_data$group<-ifelse(substr(rownames(survival_gene_data),14,15)<10,'Tumor','Normal')
survival_gene_data_tumor<-survival_gene_data[survival_gene_data$group=='Tumor',]
#universal COX analysis
library(ggplot2)
library(survival)
library(survminer)

for (i in 1:10){
  survival_gene_data_tumor$group<-ifelse(survival_gene_data_tumor[,i]>median(survival_gene_data_tumor[,i]),'high','low')
  fit <- survfit(Surv(as.numeric(OS.time), as.numeric(OS)) ~survival_gene_data_tumor$group, data = survival_gene_data_tumor)
  
  p=ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               legend.labs = c("High", "Low"), 
               legend.title = "Group", 
               palette = c( "#ff0000","#0000ff"))
  
  print(p$plot)
  ggsave(filename = paste0('K_M/',colnames(survival_gene_data_tumor)[i],'.png'))

}

#construction of prognostic model using driver_gene
source("doCOX.R")
uni_cox_df<-processUniCOX(driver_genes,survival_gene_data)

#7.lasso回归
coefficient<-dolasso(uni_cox_df,survival_gene_data)
active.Index<-which(as.numeric(coefficient)!=0)
active.coefficients<-as.numeric(coefficient)[active.Index]
sig_gene_multi_cox<-rownames(coefficient)[active.Index]

#8. 建立多因素COX回归
#perform the multi-variates cox regression using qualified genes
multi_variate_cox_2<-processMultiCOX(sig_gene_multi_cox, survival_gene_data)
ph_hypo_table=doPH_hypothesis(sig_gene_multi_cox, survival_gene_data)

#9.check the co_linearity between samples.

corration<-cor(survival_gene_data[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]],method='pearson')
library(GGally)
ggpairs(survival_gene_data[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]],
        axisLabels = 'show')+
  theme_bw()+
  theme(
    panel.background = element_rect(colour = 'black',size=0.3,fill = 'white'),
    panel.grid = element_blank()
  )
ggsave(filename = 'cor_cox.png',width =5,height = 3)

library(rms)    
vif<-rms::vif(multi_variate_cox_2)

#some people said if the square root of VIF>2, the might be co_liner.
sqrt(vif)<2
#forest plot
ggforest(model = multi_variate_cox_2,data = survival_gene_data,main = 'Hazard ratios of candidate gene',fontsize = 1)
ggsave(filename = 'multi_cox_forestplot.png',width = 7,height = 3)

#9. cox模型性能评价
#（1）C-index
C_index<-multi_variate_cox_2$concordance['concordance']

#(2) time-dependent ROC curve
#for riskscore
candidate_genes_for_cox2<-c(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05])
save(candidate_genes_for_cox2,multi_variate_cox_2,survival_gene_data,file = 'for_validatioin.Rdata')

risk_score_table_multi_cox2<-riskscore(survival_gene_data,candidate_genes_for_cox2,multi_variate_cox_2)
risk_score_table_multi_cox2<-risk_score_table_multi_cox2[risk_score_table_multi_cox2$OS.time>0,]

#timeROC
library(timeROC)
ROC.DSST<-timeROC(T=risk_score_table_multi_cox2$OS.time,delta=risk_score_table_multi_cox2$OS,
                  marker=risk_score_table_multi_cox2$total_risk_score,cause=1,
                  weighting="cox",
                  times=c(365*1,365*3,365*5,365*7),ROC=TRUE)


plot(ROC.DSST,time=1*365)        
plot(ROC.DSST,time=3*365,add=TRUE,col="blue") 
plot(ROC.DSST,time=5*365,add=TRUE,col="grey50") 

#evaluate AUCs between 3-5 years.

for_multi_ROC<-multi_ROC(time_vector = c(365*seq(3,5,1)),risk_score_table = risk_score_table_multi_cox2)

#visualization for the ROC curves of multiple time points
#plot ROC
tiff(file='AUC2.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
year3ROC<-for_multi_ROC[for_multi_ROC$Time_point==3*365,]
year4ROC<-for_multi_ROC[for_multi_ROC$Time_point==4*365,]
year5ROC<-for_multi_ROC[for_multi_ROC$Time_point==5*365,]
TP_3year = year3ROC$True_positive
FP_3year = year3ROC$False_positive
TP_4year = year4ROC$True_positive
FP_4year = year4ROC$False_positive
TP_5year = year5ROC$True_positive
FP_5year = year5ROC$False_positive
time_ROC_df=data.frame(TP_3year,FP_3year,TP_4year,FP_4year,TP_5year,FP_5year)
library(ggplot2)
p=ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", year3ROC$AUC)), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 4 years = ", sprintf("%.3f", year4ROC$AUC)), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", year5ROC$AUC)), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

print(p)
dev.off()

#K-M analysis
AUC_max<-max(for_multi_ROC$AUC)
AUC_max_time<-for_multi_ROC$Time_point[which(for_multi_ROC$AUC==AUC_max)]
AUC_max_time<-AUC_max_time[!duplicated(AUC_max_time)]
AUC_max_time<-AUC_max_time[length(AUC_max_time)]
for_multi_ROC$Time_point<-as.factor(for_multi_ROC$Time_point)
optimal_time_ROC_df<-for_multi_ROC[which(for_multi_ROC$Time_point==AUC_max_time),]
cut.off<-optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive - optimal_time_ROC_df$False_positive)]
#cut.off<-median(risk_score_table_multi_cox2$total_risk_score)
high_low<-risk_score_table_multi_cox2$total_risk_score>cut.off
high_low[high_low==TRUE]<-'high'
high_low[high_low==FALSE]<-'low'
risk_score_table_multi_forKM<-cbind(risk_score_table_multi_cox2,high_low)

#KM_plot 
library(survminer)
risk_score_table_multi_forKM$OS[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-0
risk_score_table_multi_forKM$OS.time[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-AUC_max_time
fit_km<-survfit(Surv(OS.time,OS)~high_low,data = risk_score_table_multi_forKM)
tiff(file='tmethyscore2.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
p=ggsurvplot(fit_km,conf.int = F,pval = T,legend.title="Total risk score",
             legend.labs=c('high','low'),risk.table=T,
             #legend.labs=c(paste0('>',as.character(round(cut.off,2))),
             # paste0('<=',as.character(round(cut.off,2)))),risk.table=T,
             palette=c('red','blue'),surv.median.line='hv')
print(p)
dev.off()

#predictive for patients characteristics
TCGA_pheno_data_filt2$sampleid<-substr(TCGA_pheno_data_filt2$sampleid,1,15)
rownames(TCGA_pheno_data_filt2)<-TCGA_pheno_data_filt2$sampleid
use_sample<-intersect(TCGA_pheno_data_filt2$sampleid,risk_score_table_multi_cox2$sample)
TCGA_pheno_data_filt2$sample<-TCGA_pheno_data_filt2$sampleid
TCGA_pheno_data_filt2$sampleid<-NULL

clinical_survival <- merge(risk_score_table_multi_cox2, TCGA_pheno_data_filt2, by = "sample")
data<-clinical_survival

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
#for T stage
data<-data[data$Group=="Tumor",]
T_data<-data.frame(data$T,data$total_risk_score) 
colnames(T_data)<-c('T stage','level of risk score')

rownames(T_data)<-data$sample
T_data$`T stage`<-as.factor(T_data$`T stage`)
mydata<-T_data %>% gather(key = "T stage",value = "level of risk score") %>% na.omit() 
my_comparisons <- list( c("T1", "T2"),c("T1","T3"),c("T1","T4"),c("T2","T3"),c("T2","T4"),c("T3","T4"))

tiff(file='T_stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
p<-ggboxplot(mydata,x='T stage',y='level of risk score',color='T stage',palette = 'jama',add = 'jitter')
p+stat_compare_means(comparisons = my_comparisons)
dev.off()


#for M stage
M_data<-data.frame(data$M,data$total_risk_score) #让分数变为正数
colnames(M_data)<-c('M stage','level of risk score')
rownames(M_data)<-data$sample
M_data$`M stage`<-as.factor(M_data$`M stage`)
mydata<-M_data %>% gather(key = "M stage",value = "level of risk score") %>% na.omit() 
my_comparisons <- list( c("M0", "M1"),c("M0", "MX"),c("M1", "MX"))
tiff(file='M_stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
p<-ggboxplot(mydata,x='M stage',y='level of risk score',color='M stage',palette = 'jama',add = 'jitter')
p+stat_compare_means(comparisons = my_comparisons)
dev.off()


#for N stage
N_data<-data.frame(data$N,data$total_risk_score+max(abs(data$total_risk_score))) #让分数变为正数
colnames(N_data)<-c('N stage','level of risk score')
rownames(N_data)<-data$sample
N_data$`N stage`<-as.factor(N_data$`N stage`)
mydata<-N_data %>% gather(key = "N stage",value = "level of risk score") %>% na.omit() 
my_comparisons <- list( c("N0", "N1"),c("N0", "N2"),c("N1", "N2"))
tiff(file='N_stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
p<-ggboxplot(mydata,x='N stage',y='level of risk score',color='N stage',palette = 'jama',add = 'jitter')
p+stat_compare_means(comparisons = my_comparisons)
dev.off()

#for stage
stage_data<-data.frame(data$Stage,data$total_risk_score) #让分数变为正数
colnames(stage_data)<-c('stage','level of risk score')
rownames(stage_data)<-data$sample
stage_data$`stage`<-as.factor(stage_data$`stage`)
mydata<-stage_data %>% gather(key = "stage",value = "level of risk score") %>% na.omit() 
my_comparisons <- list( c("I","II"), c("I","III"), c("II","III"), c("II","IV"), c("III","IV"))
tiff(file='stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
p<-ggboxplot(mydata,x='stage',y='level of risk score',color='stage',palette = 'jama',add = 'jitter')
p+stat_compare_means(comparisons = my_comparisons)
dev.off()

