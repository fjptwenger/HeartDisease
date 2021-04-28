
processMethylMix <- function(MET_cancer, GE_cancer, MET_normal) {
	#BiocManager::install('MethylMix')
	library(MethylMix)
  library(ggplot2)
	library(doParallel)
	cl<-makeCluster(5)
	registerDoParallel(cl)	
  
	MET_result<-MethylMix(MET_cancer_new, GE_cancer_new, MET_normal_new)

	driver_genes<-MET_result$MethylationDrivers
	#methylation-driver gene
	if(!file.exists("./cancergeneplot/")){
		dir.create("./cancergeneplot/")
	}
	for (i in driver_genes){
	  p= MethylMix_PlotModel(i,MET_result, MET_cancer, METnormal = MET_normal)
	  ggsave(filename = paste0('cancergeneplot/',i,'.png'),height = 4,width = 8)
	}

	#correlation of driver gene
	if(!file.exists("./correlationplot/")){
		dir.create("./correlationplot/")
	}
	for (i in driver_genes){
	  p= MethylMix_PlotModel(i, MET_result, MET_cancer, GE_cancer, MET_normal)
	  ggsave(filename = paste0('correlationplot/',i,'.png'),height = 4,width = 8)
	}

	#自画correlation of driver gene
	
	if(!file.exists("./correlationplot_pvalue/")){
		dir.create("./correlationplot_pvalue/")
	}
	
	for (genename in driver_genes){
	  x<-MET_cancer[genename,]
	  y<-GE_cancer[genename,]
	  corT <- cor.test(x,y)
	  z <- lm(y~x)
	  estimate <- corT$estimate
	  cor <- round(estimate,3)
	  pvalue <- corT$p.value   #以小数的方式保留位数
	  pval <- signif(pvalue,4) #以科学计数法的方式保留位数
	  pval <- format(pval,scientific=T)
	  picname=paste(genename,"cor.tiff",sep='-')
	  outTiff <- paste('correlationplot_pvalue',picname,sep='/')
	  tiff(file=outTiff,width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
	  plot(x,y,type="p",pch=16,main=paste("Cor=",cor,"(p-value=",pval,")",sep=""),xlab=paste(genename,"methylation"),ylab=paste(genename,"expression"))
	  lines(x,fitted(z),col=2)
	  dev.off()
	  line=paste(genename,round(corT$estimate,4),corT$p.value,sep='\t')
	 
	  write.table(line, "cor_result.xls",col.names = F, append = T,row.names = F,quote = F)
	}
	
	 return(MET_result)
}