library(iMKT)
library(data.table)
mktByGene <- function(data=NULL,geneList=NULL,test=NULL,population=NULL,cutoff=0.025){

	tmp <- list()
	output <- list()

	for(iter in 1:length(geneList)){
	# for(iter in 1:length(genes)){
		## Subset data by geneda
		print(iter)
		subsetGene <- data[data[['id']] == geneList[iter],]
		# subsetGene <- data[data[['id']] == genes[iter],]
		if(dim(subsetGene)[1] > 1){
			subsetGene <- subsetGene[1]
		}

		if(subsetGene$m0 == 0 | subsetGene$mi == 0 | subsetGene$p0 == 0 | subsetGene$pi == 0 | subsetGene$di == 0 | subsetGene$d0 == 0){
			tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test,stringsAsFactors=F)
			tmp[[iter]] <- tmpDf
		}
		else{
			## Set counters to 0
			# pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
			# p0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
			# mi <- 0; m0 <- 0
			# di <- 0; d0 <- 0

			## Group genes
			f <- seq(0.025,0.975,0.05)
			pi <- do.call(rbind,lapply(as.character(subsetGene[['daf0f']]), function(z) unlist(strsplit(z, split = ";")) %>% as.numeric())) %>% colSums()
			p0 <- do.call(rbind,lapply(as.character(subsetGene[['daf4f']]), function(z) unlist(strsplit(z, split = ";")) %>% as.numeric())) %>% colSums()

			mi <- subsetGene[['mi']]
			m0 <- subsetGene[['m0']]
			di <- subsetGene[['di']]
			d0 <- subsetGene[['d0']]

			## Proper formats
			daf <- cbind(f, pi, p0) %>% as.data.frame()
			names(daf) <- c("daf","Pi","P0")
			div <- cbind(mi, di, m0, d0) %>% as.data.frame()
			names(div) <- c("mi","Di","m0","D0")
			
			if(test == 'standardMKT'){
				mkt <- standardMKT(daf=daf,div=div)
				alpha <- mkt$alpha$alpha
				pvalue <- mkt$alpha$Fisher.Test
				tmpDf <- data.frame('id'=subsetGene[['id']],'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test,stringsAsFactors=F)
				# r <- c(subsetGene[['id']],population,alpha,pvalue,test)
				tmp[[iter]] <- tmpDf
				# tmp <- rbind(tmp,tmpDf)
			}else if(test == 'FWW'){

				checkDaf <- daf[daf[['daf']] > cutoff, ]

				if(sum(checkDaf[['Pi']]) == 0 | sum(checkDaf[['P0']]) == 0){
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=paste0(test))
					tmp[[iter]] <- tmpDf
					# tmp <- rbind(tmp,tmpDf)
				}else{
					fayDaf<-daf
					fayDaf[1,]<-c(0.025,0,0)
					fayDaf[2,]<-c(0.075,0,0)

					mkt <- FWW(daf=fayDaf,div=div,listCutoffs=cutoff,plot=FALSE)
					alpha <- mkt$alphaCorrected$alpha
					pvalue <- mkt$alphaCorrected$pvalue
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test,stringsAsFactors=F)
					tmp[[iter]] <- tmpDf
					# tmp <- rbind(tmp,tmpDf)
				}
			}else if(test == 'eMKT'){

				PiMinus     = sum(daf[daf[['daf']] <= cutoff,'Pi'])
				PiGreater   = sum(daf[daf[['daf']] > cutoff,'Pi'])
				P0Minus     = sum(daf[daf[['daf']] <= cutoff,'P0'])
				P0Greater   = sum(daf[daf[['daf']] > cutoff,'P0'])

				if(PiGreater== 0 | P0Greater  == 0){
				
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test,stringsAsFactors=F)
					tmp[[iter]] <- tmpDf
					
				}else{

	
					ratioP0     = P0Minus/P0Greater
					deleterious = PiMinus - (PiGreater * ratioP0)

					alpha = 1 - (((sum(daf$Pi) - deleterious)/sum(daf$P0))*(d0/di))
					m      = matrix(c(sum(daf$P0),(sum(daf$Pi) - deleterious),div[['D0']],div[['Di']]), ncol=2)
					pvalue = fisher.test(round(m))$p.value


					# mkt <- eMKT(daf=daf,div=div,cutoff=cutoff,plot=FALSE)
					# alpha <- mkt$alphaCorrected$alphaCorrected2
					# pvalue <- mkt$alphaCorrected$fisherTest
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test,stringsAsFactors=F)
					tmp[[iter]] <- tmpDf
					# tmp <- rbind(tmp,tmpDf)

				}
			}else if(test == 'aMKT'){

				mkt <- tryCatch({aMKT(daf=daf,div=div,xlow=0.1,xhigh=0.9,plot=F)},error=function(e){mkt<-NULL})

				if(is.null(mkt)){
					alpha <- NA
					pvalue <- NA
				}else{
					alpha <- mkt$alpha <- mkt$alphaCorrected$alphaAsymptotic
					if(alpha > 0 & alpha < 1  & mkt$alphaCorrected$ciLow > 0 & !is.nan(mkt$alphaAsymptotic$CI_low)){
						pvalue <- 0.005
					}else if(alpha < 0 & mkt$alphaCorrected$ciHigh < 0 & !is.nan(mkt$alphaAsymptotic$ciHigh)){
						pvalue <- 0.005
					}else{
						pvalue <- 1
					}
				}

				tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test,stringsAsFactors=F)
				tmp[[iter]] <- tmpDf
				# tmp <- rbind(tmp,tmpDf)
			}else if(test == 'caMKT'){

				mkt <- tryCatch({aMKT(daf=daf,div=div,xlow=0,xhigh=1,plot=FALSE)},error=function(e){mkt<-NULL})

				if(is.null(mkt)){
					alpha <- NA
					pvalue <- NA
				}else{
					alpha <- mkt$alpha <- mkt$alphaCorrected$alphaAsymptotic
					if(alpha > 0 & alpha < 1  & mkt$alphaCorrected$ciLow > 0 & !is.nan(mkt$alphaAsymptotic$CI_low)){
						pvalue <- 0.005
					}else if(alpha < 0 & mkt$alphaCorrected$ciHigh < 0 & !is.nan(mkt$alphaAsymptotic$ciHigh)){
						pvalue <- 0.005
					}else{
						pvalue <- 1
					}
				}


				tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test,stringsAsFactors=F)
				tmp[[iter]] <- tmpDf
				# tmp <- rbind(tmp,tmpDf)
			}
		}
	}
	## Delete list names to not overwrite current content (saving list positions as populations each gene)
	tmp <- rbindlist(tmp)
	# tmp[['pvalue']] <- p.adjust(tmp[['pvalue']],method='fdr')

	analyzable <- matrix(c(dim(tmp %>% na.omit)[1],mean(tmp[['alpha']],na.rm=T),sd(tmp[['alpha']],na.rm=T)),ncol=3,nrow=1)
	positive <- matrix(c(dim(tmp[alpha > 0 & pvalue < 0.05])[1],mean(tmp[alpha > 0 & pvalue < 0.05,alpha],na.rm=T),sd(tmp[alpha>0 & pvalue < 0.05,alpha],na.rm=T)),ncol=3,nrow=1)
	negative <- matrix(c(dim(tmp[alpha<0 & pvalue < 0.05])[1],mean(tmp[alpha<0 & pvalue < 0.05,alpha],na.rm=T),sd(tmp[alpha<0 & pvalue < 0.05,alpha],na.rm=T)),ncol=3,nrow=1)

	# Formating alpha table with total, mean and sd
	output[['alphaTable']] <- as.data.frame(rbind(analyzable,positive,negative))
	colnames(output[['alphaTable']]) <- c('N','mean','sd')
	output[['alphaTable']][['test']] <- test
	output[['alphaTable']][['type']] <- c('analyzable','positive','negative')
	output[['alphaTable']][['pop']] <- population

	output[['alphaTable']] <- output[['alphaTable']][,c('test','type','N','mean','sd','pop')]
	
	output[['posSigGenes']] <- as.character(tmp[alpha>0 & pvalue < 0.05,id])

	return(output)
}


dosByGene <- function(data=NULL,geneList=NULL,population=NULL,cutoff=NULL){

	tmp <- data.frame('id'=character(),'population'=character(),'dos'=integer(),'pvalue'=integer())
	output <- list()
	for(iter in 1:length(geneList)){
		## Subset data by geneda
		print(iter)
		subsetGene <- data[data[['id']] == geneList[iter],]

		if(dim(subsetGene)[1] > 1){
			subsetGene <- subsetGene[1]
		}

		if(subsetGene$m0 == 0 | subsetGene$mi == 0 | subsetGene$p0 == 0 | subsetGene$pi == 0 | subsetGene$di == 0 | subsetGene$d0 == 0){
			tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'dos'=NA,'pvalue'=NA)
			tmp <- rbind(tmp,tmpDf)
		}
		else{
			# DoS = D(n)/(D(n) + D(s)) - P(n)/(P(n) + P(s))
			if(is.null(cutoff)){
				dos <- (subsetGene$di/(subsetGene$di + subsetGene$d0)) - (subsetGene$pi/(subsetGene$pi + subsetGene$p0))

				m <- matrix(c(subsetGene$p0, subsetGene$pi, subsetGene$d0, subsetGene$di),ncol = 2)
				pvalue <- fisher.test(m)$p.value
				
				tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'dos'=dos,'pvalue'=pvalue)
				tmp <- rbind(tmp,tmpDf)

			}
			else{

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
				## Cleaning slightly deleterious mutation by cutoff
				Pi <- sum(daf$Pi); P0 <- sum(daf$P0)
				dafBelowCutoff <- daf[daf$daf <= cutoff, ]
				dafAboveCutoff <- daf[daf$daf > cutoff, ]


				mktTableStandard <- data.frame(Polymorphism = c(sum(daf$P0),sum(daf$Pi)), divergence = c(div$D0, div$Di),row.names = c("Neutral class", "Selected class"))
											   
				mktTable <- data.frame(`dafBelowCutoff` = c(sum(dafBelowCutoff$P0),sum(dafBelowCutoff$Pi)), `dafAboveCutoff` = c(sum(dafAboveCutoff$P0),
					sum(dafAboveCutoff$Pi)), row.names = c("neutralClass", "selectedClass"))

				fNeutral <- mktTable[1, 1]/sum(daf$P0)
				piNeutralBelowCutoff <- Pi * fNeutral
				piWd <- mktTable[2, 1] - piNeutralBelowCutoff
				piNeutral <- round(piNeutralBelowCutoff + mktTable[2,2])

				if(piNeutral == 0){
					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'dos'=NA,'pvalue'=pvalue)
					tmp <- rbind(tmp,tmpDf)
				}
				else{
					dos <- (subsetGene$d0/(subsetGene$d0 + subsetGene$di)) - (subsetGene$p0/(subsetGene$p0 + piNeutral))
					m <- matrix(c(subsetGene$p0, subsetGene$pi, subsetGene$d0, subsetGene$di),ncol = 2)
					pvalue <- fisher.test(m)$p.value

					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'dos'=dos,'pvalue'=pvalue)
					tmp <- rbind(tmp,tmpDf)
				}
			}
		}
	}
	## Delete list names to not overwrite current content (saving list positions as populations each gene)
	tmp <- as.data.table(tmp)
	# tmp[['pvalue']] <- p.adjust(tmp[['pvalue']],method='fdr')
	output[['dos']] <- tmp

	allGenes <- matrix(c(dim(tmp %>% na.omit)[1],mean(tmp[['dos']],na.rm=T),sd(tmp[['dos']],na.rm=T)),ncol=3,nrow=1)
	positive <- matrix(c(dim(tmp[dos>0 & pvalue < 0.05])[1],mean(tmp[dos>0 & pvalue < 0.05,dos],na.rm=T),sd(tmp[dos>0 & pvalue < 0.05,dos],na.rm=T)),ncol=3,nrow=1)
	negative <- matrix(c(dim(tmp[dos<0 & pvalue < 0.05])[1],mean(tmp[dos<0 & pvalue < 0.05,dos],na.rm=T),sd(tmp[dos<0 & pvalue < 0.05,dos],na.rm=T)),ncol=3,nrow=1)

	# Formating alpha table with total, mean and sd
	output[['dosTable']] <- as.data.frame(rbind(allGenes,positive,negative))

	colnames(output[['dosTable']]) <- c('N','mean','sd')

	if(is.null(cutoff)){
		output[['dosTable']][['Pi']] <- 'rawPi'
		output[['dosTable']][['type']] <- c('allGenes','positive','negative')
		output[['dosTable']][['pop']] <- population

		output[['dosTable']] <- output[['dosTable']][,c('Pi','type','N','mean','sd','pop')]
		
		output[['posSigGenes']] <- as.character(tmp[dos>0 & pvalue < 0.05,id])
	}else{
		output[['dosTable']][['Pi']] <- 'PiNeutral'
		output[['dosTable']][['type']] <- c('allGenes','positive','negative')
		output[['dosTable']][['pop']] <- population

		output[['dosTable']] <- output[['dosTable']][,c('Pi','type','N','mean','sd','pop')]
		
		output[['posSigGenes']] <- as.character(tmp[dos>0 & pvalue < 0.05,id])
	}

	return(output)
}


# else if(test == 'FWW' & cutoff==0.1){

# 				checkDaf <- daf[daf[['daf']] > cutoff, ]

# 				if(sum(checkDaf[['Pi']]) == 0 | sum(checkDaf[['P0']]) == 0){
# 					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test,stringsAsFactors=F)
# 					tmp <- rbind(tmp,tmpDf)
# 				}else{
# 					daf[1,]<-c(0.025,0,0)
# 					daf[2,]<-c(0.075,0,0)
# 					mkt <- FWW(daf=daf,div=div,listCutoffs=cutoff,plot=FALSE)
# 					alpha <- mkt$alphaCorrected$alpha
# 					pvalue <- mkt$alphaCorrected$pvalue
# 					tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test,stringsAsFactors=F)
# 					tmp[[iter]] <- tmpDf
# 					tmp <- rbind(tmp,tmpDf)
# 				}
# 			}else if(test == 'eMKT' & cutoff==0.1){

			# 	if(sum(daf[['Pi']]) == 0 | sum(daf[['P0']]) == 0){
			# 		tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=NA,'pvalue'=NA,'test'=test,stringsAsFactors=F)
			# 		tmp <- rbind(tmp,tmpDf)
			# 	}else{

			# 		mkt <- eMKT(daf=daf,div=div,cutoff=cutoff,plot=FALSE)
			# 		alpha <- mkt$alphaCorrected$alphaCorrected2
			# 		pvalue <- mkt$alphaCorrected$fisherTest
			# 		tmpDf <- data.frame('id'=subsetGene$id,'pop'=population,'alpha'=alpha,'pvalue'=pvalue,'test'=test,stringsAsFactors=F)
			# 		tmp[[iter]] <- tmpDf
			# 		# tmp <- rbind(tmp,tmpDf)
			# 	}
			# }