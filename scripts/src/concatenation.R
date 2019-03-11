library(splitstackshape)

samplingGenes <- function(geneList,B,bins,seed=213,path=NULL){
	set.seed(seed)
	
	bootSamples <- matrix(sample(geneList,size = B*bins,replace=T),B,bins)
	
	output <- list()
	for (i in 1:B) {
		output[[i]] <- bootSamples[i,]
	}

	if(!is.null(path)){
		save(output,file=paste0(path,'bootSampleBin',b,'.RData'))
	}

	return(output)
}


sampleAnalysis <- function(data,sampling,bins,population){
	data <- data %>% as.data.table
	data <- data[Pop==population]
	data$Name <- as.character(data$Name)
	output <- data.frame()

	for (i in 1:length(sampling)) {
		print(i)	
		sampling[[i]] <- sampling[[i]] %>% as.data.frame
		colnames(sampling[[i]]) <- 'Name'
		# subsetData <- data[match(data$Name, sampling[[i]]),] %>% na.omit()
		subsetData <- merge(sampling[[i]],data,by='Name',all.x=T)

		divergence <- subsetData[,c('mi', 'di', 'm0', 'd0')]
		colnames(divergence) <- c('mi', 'Di', 'm0', 'D0')
		divergence <- divergence %>% summarize(mi=sum(mi),Di=sum(Di),m0=sum(m0),D0=sum(D0))

		sfs0f <- as.data.frame(subsetData[['DAF0f']])
		colnames(sfs0f) <- 'sfs0f'
		sfs0f <- cSplit(sfs0f, 'sfs0f', ';') %>% colSums
		sfs0f <- as.numeric(sfs0f)

		sfs4f <- as.data.frame(subsetData[['DAF4f']]) 
		colnames(sfs4f) <- 'sfs4f'
		sfs4f <- cSplit(sfs4f, 'sfs4f', ';') %>% colSums
		sfs4f <- as.numeric(sfs4f)

		daf <- cbind(daf=c(0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975), sfs0f, sfs4f) %>% as.data.frame
		colnames(daf) <- c('daf','Pi','P0')
		# daf<-daf[daf$daf <= 0.9, ]
		D0 <- divergence$D0 %>% sum
		Di <- divergence$Di %>% sum
		Pi <- sum(daf$Pi)
		P0 <- sum(daf$P0)

		if (P0 == 0 | D0  == 0) {
			tmp <- data.frame(NA,NA,NA,NA,NA,NA,NA)
			output <- rbind(output,tmp)
		}else{
			resultStandard <- standardMKT(daf,divergence)
			alphaStandard <- resultStandard$alpha.symbol

			resultDGRP0.05 <- DGRP(daf,divergence,listCutoffs = c(0.025))
			alphaDGRP0.05 <- resultDGRP0.05$Results$alpha.symbol

			resultDGRP0.15 <- DGRP(daf,divergence,listCutoffs = c(0.075))
			alphaDGRP0.15 <- resultDGRP0.15$Results$alpha.symbol

			resultFWW0.05 <- FWW(daf,divergence,listCutoffs = c(0.025))
			alphaFWW0.05 <- resultFWW0.05$Results$alpha.symbol

			resultFWW0.15 <- FWW(daf,divergence,listCutoffs = c(0.075))
			alphaFWW0.15 <- resultFWW0.15$Results$alpha.symbol

			resultiMK0.1 <- NULL
			alphaAsymptotic0.1 <- NULL
			resultiMK <- NULL
			alphaAsymptotic <- NULL

		}

		tmp <- data.frame(alphaStandard,alphaDGRP0.05,alphaDGRP0.15,alphaFWW0.05,alphaFWW0.15)
		# tmp <- data.frame(alphaStandard,alphaDGRP0.05,alphaDGRP0.15,alphaFWW0.05,alphaFWW0.15,alphaAsymptotic,alphaAsymptotic0.1)

		output <- rbind(output,tmp)	
	}
	output[['bin']] <- bins
	return(output)
}