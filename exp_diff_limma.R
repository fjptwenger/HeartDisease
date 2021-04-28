
processData <- function(exprset) {
	suppressMessages(library(edgeR))
	suppressMessages(library(limma))
	exprset<-COAD_READ_expr
	group_list=ifelse(as.numeric(substr(colnames(exprset),14,15))<10,'Tumor','Normal')
	sample_type<-factor(group_list)	
	
	design<-model.matrix(~0+sample_type)
	colnames(design)<-levels(sample_type)
	rownames(design)<-colnames(exprset)
	dge <- DGEList(counts=exprset)
	dge<-calcNormFactors(dge)	
	v<-voom(dge,design,normalize='quantile',)
	fit<-lmFit(v,design)
	constrasts<-paste(rev(levels(sample_type)),collapse='-')	
	cont.matrix<-makeContrasts(contrast=constrasts,levels=sample_type)
	fit2<-contrasts.fit(fit,cont.matrix)
	fit2<-eBayes(fit2)
	suppressMessages(library(dplyr))
	suppressMessages(library(ggplot2))
	DEG=topTable(fit2,coef=constrasts, n=Inf) %>% na.omit()
		
	DEG$change=as.factor(
		ifelse(
			DEG$P.Value <0.05 & abs(DEG$logFC)>0.585,
				ifelse(DEG$logFC >0.585,'UP','DOWN'),'NOT'
		)
	)	
	return(DEG)
}


