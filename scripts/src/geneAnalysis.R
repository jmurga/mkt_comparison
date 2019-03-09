mktByGene <- function(data=NULL,geneList=NULL,test=NULL,population=NULL,cutoff=NULL){

	tmp <- data.frame('id'=character(),'population'=character(),'alpha'=integer(),'pvalue'=integer(),test=character())
	output <- list()
	for(iter in 1:length(geneList)){
		## Subset data by geneda
		print(iter)
		subsetGene <- data[data[['id']] == geneList[iter],]
		if(dim(subsetGene)[1] > 1){
			subsetGene <- subsetGene[1]
		}

		if(subsetGene$m0 == 0 | subsetGene$mi == 0 | subsetGene$p0 == 0 | subsetGene$pi == 0 | subsetGene$di == 0 | subsetGene$d0 == 0){
			tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test)
			tmp <- rbind(tmp,tmpDf)
		}
		else{
			## Set counters to 0
			pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
			p0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
			f <- seq(0.025,0.975,0.05)
			mi <- 0; m0 <- 0
			di <- 0; d0 <- 0

			## Group genes
			pi <- do.call(rbind,lapply(as.character(subsetGene$DAF0f), function(z) unlist(strsplit(z, split = ";")) %>% as.numeric())) %>% colSums()
			p0 <- do.call(rbind,lapply(as.character(subsetGene$DAF4f), function(z) unlist(strsplit(z, split = ";")) %>% as.numeric())) %>% colSums()
			mi <- subsetGene$mi
			m0 <- subsetGene$m0
			di <- subsetGene$di
			d0 <- subsetGene$d0

			## Proper formats
			daf <- cbind(f, pi, p0); daf <- as.data.frame(daf)
			names(daf) <- c("daf","Pi","P0")
			div <- cbind(mi, di, m0, d0); div <- as.data.frame(div)
			names(div) <- c("mi","Di","m0","D0")
			
			if(test == 'standardMKT'){
				mkt <- standardMKT(daf=daf,div=div)
				alpha <- mkt$alpha.symbol
				pvalue <- mkt$`Fishers exact test P-value
				tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test)
				tmp <- rbind(tmp,tmpDf)				`
			}
			else if(test == 'FWW' & cutoff==0.05){
				checkDaf <- daf[daf[['daf']] > cutoff, ]

				if(sum(checkDaf[['Pi']]) == 0 | sum(checkDaf[['P0']]) == 0){
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}else{
					daf[1,]<-c(0.025,0,0)
					mkt <- FWW(daf=daf,div=div,listCutoffs=cutoff,plot=FALSE)
					alpha <- mkt$Results$alpha.symbol
					pvalue <- mkt$Results$`Fishers exact test P-value`
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}
			}
			else if(test == 'FWW' & cutoff==0.1){

				checkDaf <- daf[daf[['daf']] > cutoff, ]

				if(sum(checkDaf[['Pi']]) == 0 | sum(checkDaf[['P0']]) == 0){
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}
				else{
					daf[1,]<-c(0.025,0,0)
					daf[2,]<-c(0.075,0,0)
					mkt <- FWW(daf=daf,div=div,listCutoffs=cutoff,plot=FALSE)
					alpha <- mkt$Results$alpha.symbol
					pvalue <- mkt$Results$`Fishers exact test P-value`	
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}
			}
			else if(test == 'eMKT' & cutoff==0.05){

				## Cleaning slightly deleterious mutation by cutoff
				pi <- sum(pi); p0 <- sum(p0)
				dafBelowCutoff <- daf[daf$daf <= cutoff, ]
				dafAboveCutoff <- daf[daf$daf > cutoff, ]

				mktTableStandard <- data.frame(Polymorphism = c(sum(daf$p0), 
					sum(daf$pi)), divergence = c(div$d0, div$di), 
					row.names = c("Neutral class", "Selected class"))
											   
				mktTable <- data.frame(`dafBelowCutoff` = c(sum(dafBelowCutoff$p0), 
					sum(dafBelowCutoff$pi)), `dafAboveCutoff` = c(sum(dafAboveCutoff$p0), 
					sum(dafAboveCutoff$pi)), row.names = c("neutralClass", 
					"selectedClass"))
											   
				fNeutral <- mktTable[1, 1]/sum(daf$p0)
				piNeutralBelowCutoff <- pi * fNeutral
				piWd <- mktTable[2, 1] - piNeutralBelowCutoff
				piNeutral <- round(piNeutralBelowCutoff + mktTable[2, 
					2])

				if(sum(dafAboveCutoff[['Pi']]) == 0 | piNeutral == 0){
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}
				else{
					
					alpha <- 1 - ((piNeutral/p0) * (mktTableStandard[1, 2]/mktTableStandard[2, 2]))
					
					m <- matrix(c(p0, piNeutral, div$d0, div$di),ncol = 2)
					pvalue <- fisher.test(m)$p.value

					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}
			}
			else if(test == 'eMKT' & cutoff==0.1){
				checkDaf <- daf[daf[['daf']] > cutoff, ]

				if(sum(checkDaf[['Pi']]) == 0 | sum(checkDaf[['P0']]) == 0){
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}
				else{
					daf[1,]<-c(0.025,0,0)
					daf[2,]<-c(0.075,0,0)					
					mkt <- DGRP(daf=daf,div=div,listCutoffs=cutoff,plot=FALSE)
					alpha <- mkt$Results$alpha.symbol
					pvalue <- mkt$Results$`Fishers exact test P-value`	
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test)
					tmp <- rbind(tmp,tmpDf)
				}
			}
			else if(test == 'aMKT'){
				mkt <- tryCatch({asymptoticMKT(daf=daf,div=div,xlow=0.1,xhigh=0.9,plot=FALSE)},error=function(e){mkt<-NULL})
				if(is.null(mkt)){
					alpha <- NA
					pvalue <- NA
				}else{
					alpha <- mkt$Results$alpha.symbol
					pvalue <- mkt$Results$`Fishers exact test P-value`
				}

				tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test)
				tmp <- rbind(tmp,tmpDf)
			}
			else if(test == 'caMKT'){
				mkt <- tryCatch({iMKT(daf=daf,div=div,xlow=0.1,xhigh=0.9,plot=FALSE)},error=function(e){mkt<-NULL})
				if(is.null(mkt)){
					alpha <- NA
					pvalue <- NA
				}else{
					alpha <- mkt$Results$alpha.symbol
					pvalue <- mkt$Results$alpha.symbol
				}

				tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test)
				tmp <- rbind(tmp,tmpDf)
			}	
		}
	}
	## Delete list names to not overwrite current content (saving list positions as populations each gene)
	tmp <- as.data.table(tmp)
	# tmp[['pvalue']] <- p.adjust(tmp[['pvalue']],method='fdr')

	allGenes <- matrix(c(dim(tmp %>% na.omit)[1],mean(tmp[['alpha']],na.rm=T),sd(tmp[['alpha']],na.rm=T)),ncol=3,nrow=1)
	positive <- matrix(c(dim(tmp[alpha>0 & pvalue < 0.05])[1],mean(tmp[alpha>0 & pvalue < 0.05,alpha],na.rm=T),sd(tmp[alpha>0 & pvalue < 0.05,alpha],na.rm=T)),ncol=3,nrow=1)
	negative <- matrix(c(dim(tmp[alpha<0 & pvalue < 0.05])[1],mean(tmp[alpha<0 & pvalue < 0.05,alpha],na.rm=T),sd(tmp[alpha<0 & pvalue < 0.05,alpha],na.rm=T)),ncol=3,nrow=1)

	# Formating alpha table with total, mean and sd
	output[['alphaTable']] <- as.data.frame(rbind(allGenes,positive,negative))
	colnames(output[['alphaTable']]) <- c('N','mean','sd')
	output[['alphaTable']][['test']] <- test
	output[['alphaTable']][['type']] <- c('allGenes','positive','negative')
	output[['alphaTable']][['pop']] <- population

	output[['alphaTable']] <- output[['alphaTable']][,c('test','type','N','mean','sd','pop')]
	
	output[['posSigGenes']] <- as.character(tmp[alpha>0 & pvalue < 0.05,id])

	return(output)

}