

processUniCOX <- function(genelist, survival_gene_data) {
  library(survival)
  library(survminer)
  genelist<-gsub(genelist,pattern = '-',replacement = '_')
  genelist<-gsub(genelist,pattern = '/',replacement = '_OR_')
  uni_cox<-function(single_gene){
    
    formula<-as.formula(paste0('Surv(OS.time, OS) ~ ', single_gene))
    surv_uni_cox<-summary(coxph(formula,data = survival_gene_data))
    ph_hypothesis_p<-cox.zph(coxph(formula,data = survival_gene_data))$table[1,3]
    if(surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){
      single_cox_report<-data.frame(
          'uni_cox_sig_genes'=single_gene,
          #'beta'=surv_uni_cox$coefficients[,1],
          'HR'=exp(surv_uni_cox$coefficients[,1]),
          'HR.95L'=surv_uni_cox$conf.int[,3],
          'HR.95H'=surv_uni_cox$conf.int[,4], 
          'p-value'=surv_uni_cox$coefficients[,5]
          #'z_pvalue'=surv_uni_cox$coefficients[,5],
          #'wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
          #'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3])
      )
      single_cox_report
    }
    
  }
  uni_cox_list<-lapply(genelist, uni_cox)
  do.call(rbind,uni_cox_list)
}

dolasso <- function(uni_cox_df,survival_gene_data) {
	suppressMessages(library(glmnet))
	survival_gene_data<-survival_gene_data[survival_gene_data$OS.time>0,]#要删除os.time等于0的，不然会报错
	x<-as.matrix(survival_gene_data[,uni_cox_df$uni_cox_sig_genes])
	y<-survival_gene_data[,c('OS.time','OS')]
	names(y)<-c('time','status')
	y$time<-as.double(y$time)
	y$status<-as.double(y$status)
	y<-as.matrix(survival::Surv(y$time,y$status))
	lasso_fit<-cv.glmnet(x,y,family='cox',type.measure = 'deviance')
	fit<-glmnet(x,y,family = 'cox')

	coefficient<-coef(lasso_fit,s=lasso_fit$lambda.min)	
	return(coefficient)
}

doPH_hypothesis<-function(sig_gene_multi_cox, survival_gene_data){
	formula_for_multivariate<-as.formula(paste0('Surv(OS.time, OS) ~ ', paste(sig_gene_multi_cox,sep = '',collapse = '+')))
	multi_variate_cox<-coxph(formula_for_multivariate,data=survival_gene_data)
	#check if variances are supported by PH hypothesis
	ph_hypothesis_multi<-cox.zph(multi_variate_cox)
	#the last row of the table records the test results on the Global model, delete it.
	ph_hypo_table<-ph_hypothesis_multi$table[-nrow(ph_hypothesis_multi$table),]	
	return(ph_hypo_table)

}


processMultiCOX <- function(sig_gene_multi_cox, survival_gene_data) {
  library(survival)
  library(survminer)
	ph_hypo_table=doPH_hypothesis(sig_gene_multi_cox, survival_gene_data)	
	#remove variances not supported by ph hypothesis and perform the 2nd regression
	formula_for_multivariate_2<-as.formula(paste0('Surv(OS.time, OS) ~ ', paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep = '',collapse = '+')))
	multi_variate_cox_2<-coxph(formula_for_multivariate_2,data = survival_gene_data)
	return(multi_variate_cox_2)
}

riskscore<-function(survival_gene_data,candidate_genes_for_cox,cox_report){
  suppressMessages(library(dplyr))
  samplenames<-survival_gene_data$sample %>% na.omit()
  
  survival_gene_data<-survival_gene_data[survival_gene_data$sample %in% samplenames,]
  rownames(survival_gene_data)<-survival_gene_data$sample
  survival_gene_data<-survival_gene_data[samplenames,]
  
  risk_score_table<-survival_gene_data[samplenames,candidate_genes_for_cox]
  
  rownames(risk_score_table)<-survival_gene_data$sample 
  
  all_riskscore=list()
  for(each_sample in rownames(risk_score_table)){
    each_total_score=0
    for(each_sig_gene in colnames(risk_score_table)){
      
     each_total_score=each_total_score+risk_score_table[each_sample,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
    
    }
    all_riskscore=append(all_riskscore,each_total_score)
    
  }
  risk_score_table$total_risk_score<-as.numeric(all_riskscore)
  risk_score_table<-cbind(survival_gene_data[,c('sample','OS.time','OS')],risk_score_table)
  #risk_score_table<-cbind(risk_score_table,'total_risk_score'=exp(rowSums(risk_score_table))) %>%
  #  cbind(survival_gene_data[,c('title','OS.time','OS')])
  risk_score_table<-risk_score_table[,c('sample','OS.time','OS',candidate_genes_for_cox,'total_risk_score')]
 return(risk_score_table)
}

#for ROC
multi_ROC<-function(time_vector,risk_score_table){
  library(survivalROC)
  single_ROC<-function(single_time){
    for_ROC<-survivalROC(Stime = risk_score_table$OS.time,
                         status=risk_score_table$OS,
                         marker=risk_score_table$total_risk_score,
                         predict.time = single_time,method = 'KM')
    data.frame(
        'True_positive'=for_ROC$TP,'False_positive'=for_ROC$FP,
        'Cut_values'=for_ROC$cut.values,
        'Time_point'=rep(single_time,length(for_ROC$TP)),
        'AUC'=rep(for_ROC$AUC,length(for_ROC$TP))
    )
    
  }
  multi_ROC_list<-lapply(time_vector, single_ROC)
  do.call(rbind,multi_ROC_list)
}


processUniGLM <- function(genelist,drivergene_OS_mut) {
  library(rms)
   uni_glm_run<-function(single_gene){
    formula<-as.formula(paste0('OS ~ ', single_gene))
    uni_glm<-summary(glm(formula,data = drivergene_OS_mut,family = binomial(link="logit"), x=T))
    if(uni_glm$coefficients[,4][2]<0.05){
      single_report<-data.frame(
        'uni_sig_genes'=single_gene,
        'Pr(>|z|)'=as.character(uni_glm$coefficients[,4][2]),
        'z value)'=as.character(uni_glm$coefficients[,3][2]),
        'Estimate'=as.character(uni_glm$coefficients[,1][2]),
        'Std. Error'=as.character(uni_glm$coefficients[,2][2])
      )
      single_report
    }
    
  }
  uni_list<-lapply(genelist, uni_glm_run)
  do.call(rbind,uni_list)
}

processMultiglm<- function(uni_gene_glm, drivergene_OS_mut) {
  formula_for_multiglm<-as.formula(paste0('OS ~ ', paste(uni_gene_glm$uni_sig_genes,sep = '',collapse = '+')))
  multiglm_result<-glm(formula_for_multiglm,data = drivergene_OS_mut,family = binomial(link="logit"), x=T)
  return(multiglm_result)
}  

processMultilrm<- function(uni_gene_lrm, drivergene_OS_mut) {
  formula_for_multilrm<-as.formula(paste0('OS ~ ', paste(uni_gene_glm$uni_sig_genes,sep = '',collapse = '+')))
  multilrm_result<-lrm(formula_for_multilrm,data = drivergene_OS_mut,y=T, x=T)
  return(multilrm_result)
} 



