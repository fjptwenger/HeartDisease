
processMethy <- function(beta,logfc_cutoff) {

	suppressMessages(library(tibble))
	suppressMessages(library(ChAMP))

	group_list <- beta$pd$group_list
	memory.limit(2000000)#设置分配内存，不然会报错
	methy_DMP <- champ.DMP(beta = beta$beta,pheno=group_list)

	df_DMP <- methy_DMP$Tumor_to_Normal
	
	df_DMP_filt<-df_DMP[df_DMP$adj.P.Val<=0.05 & abs(df_DMP$logFC)>=logfc_cutoff,]
	df_DMP_siggene=df_DMP_filt[df_DMP_filt$gene!="",]
	df_DMP_siggene$change <-ifelse(df_DMP_siggene$logFC>logfc_cutoff,'UP','DOWN')
	return(df_DMP_siggene)	

}
