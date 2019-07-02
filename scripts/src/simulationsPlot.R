library(iMKT)
library(forcats)
library(cowplot)
library(reshape2)

mktOnsimulatedData <- function(scenario,simulationsPath){
	output <- NULL
	plotTable <- list()

	setwd(paste0(simulationsPath,scenario,'/'))
	listDaf <- lapply(list.files(pattern = 'daf',recursive=T),fread)
	listDiv <- lapply(list.files(pattern = 'div',recursive=T),fread,sep='\t',header=T)
	
	for(simulation in 1:length(listDaf)){

		# Format data from each replicate
		colnames(listDaf[[simulation]]) <- c('daf','Pi','P0')
		
		if(scenario=='length1e8'){
			listDiv[[simulation]]$m0<-1e8
			listDiv[[simulation]]$mi <- 1e8
		}
		else if(scenario=='length1e6'){
			listDiv[[simulation]]$m0<-1e6
			listDiv[[simulation]]$mi<-1e6
		}
		else{
			listDiv[[simulation]]$m0<-1e7
			listDiv[[simulation]]$mi<-1e7
		}
		
		colnames(listDiv[[simulation]]) <- c('Di', 'D0','trueAlpha','m0','mi')
		trueAlpha <- listDiv[[simulation]][,c('trueAlpha')]
		divergence <- listDiv[[simulation]][,c('Di','D0','m0','mi')]

		# Execute each mkt test through iMKT package
		## StandardMKT
		resultStandard <- standardMKT(listDaf[[simulation]],divergence)
		alphaStandard <- resultStandard$alpha.symbol

		##DGRP
		resultDGRP0.05 <- DGRP(listDaf[[simulation]],divergence,listCutoffs = c(0.025))
		alphaDGRP0.05 <- resultDGRP0.05$Results$alpha.symbol

		resultDGRP0.1 <- DGRP(listDaf[[simulation]],divergence,listCutoffs = c(0.075))
		alphaDGRP0.1 <- resultDGRP0.1$Results$alpha.symbol
		
		##FWW
		resultFWW0.05 <- FWW(listDaf[[simulation]],divergence,listCutoffs = c(0.025))
		alphaFWW0.05 <- resultFWW0.05$Results$alpha.symbol

		resultFWW0.1 <- FWW(listDaf[[simulation]],divergence,listCutoffs = c(0.075))
		alphaFWW0.1 <- resultFWW0.1$Results$alpha.symbol

		##Asymptotic
		alphaAsymptotic0.1 <- tryCatch({
			resultiMK0.1 <- iMKT(listDaf[[simulation]],divergence,xlow = 0.1, xhigh = 0.9)
			alphaAsymptotic0.1 <- resultiMK0.1$`Asymptotic MK table`$alpha_asymptotic
			},error=function(e){ alphaAsymptotic0.1 <- NA; alphaAsymptotic <- NA}
		)

		alphaAsymptotic <- tryCatch({
			resultiMK <- iMKT(listDaf[[simulation]],divergence,xlow = 0, xhigh = 1)
			alphaAsymptotic <- resultiMK$`Asymptotic MK table`$alpha_asymptotic
			},error=function(e){ alphaAsymptotic0.1 <- NA; alphaAsymptotic <- NA})

		result <- data.frame(simulation,alphaStandard,alphaDGRP0.05,alphaDGRP0.1,alphaFWW0.05,alphaFWW0.1,alphaAsymptotic,alphaAsymptotic0.1,trueAlpha)
		output <- rbind(output,result)
 
	}
	
	# Preparing data to plot
	dataPlot<-melt(output,id=c('simulation'))

	# Errors
	standardError <- mean(abs(output[['alphaStandard']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRPError <- mean(abs(output[['alphaDGRP0.05']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRP1Error <- mean(abs(output[['alphaDGRP0.1']]-output[['trueAlpha']]),na.rm=TRUE)
	FWWError <- mean(abs(output[['alphaFWW0.05']]-output[['trueAlpha']]),na.rm=TRUE)
	FWW1Error <- mean(abs(output[['alphaFWW0.1']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptoticError <- mean(abs(output[['alphaAsymptotic']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptotic1Error <- mean(abs(output[['alphaAsymptotic0.1']]-output[['trueAlpha']]),na.rm=TRUE)

	error <- matrix(c(standardError,DGRPError,DGRP1Error,FWWError,FWW1Error,asymptoticError,asymptotic1Error,0),nrow=1,ncol=8)

	meanSd <- sapply(output, function(cl) list(means=mean(cl,na.rm=TRUE), sds=sd(cl,na.rm=TRUE)))
	meanSd <- meanSd[,-1]
	meanSd <- rbind(meanSd,error)
	meanSd <- as.data.frame(meanSd)
	rownames(meanSd) <- c('mean','sd','error')

	meanSd$scenario <- scenario

	# Ploting
	plotAlpha <- ggplot(dataPlot, aes(x=variable, y=value, fill=variable)) + 
		geom_boxplot(color="grey20",alpha=0.7) + 
		labs(x = "MKT methods", y=expression(italic(α))) + 
		themePublication() + 
		scaleFillPublication(name="Method", labels=c("alphaStandard"="Standard", "alphaDGRP0.05" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW0.05" = "FWW 5%", "alphaFWW0.1" = "FWW 10%","alphaAsymptotic"="Asymptotic MKT","alphaAsymptotic0.1"="Asymptotic MKT", "trueAlpha"="Real alpha")) + 
		scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + 
		scale_x_discrete(labels=c("alphaStandard"="Standard", "alphaDGRP0.05" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW0.05" = "FWW 5%", "alphaFWW0.1" = "FWW 10%","alphaAsymptotic"="Asymptotic MKT","alphaAsymptotic0.1"="Asymptotic MKT", "trueAlpha"="Real alpha")) +  
		guides(fill=FALSE) + 
		ggtitle(paste0(scenario)) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),axis.text.y= element_text(size=20),plot.title=element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24))

	method <- rep(c("Standard","eMKT cutoff 5%","eMKT DGRP 10%","FWW cutoff 5%","MKT FWW cutoff 10%","Asymptotic MKT","Asymptotic MKT","Real alpha"),2)
					 

	plotTable[['plot']] <- plotAlpha
	plotTable[['table']] <- meanSd
	plotTable[['data']] <- dataPlot
					 
	return(plotTable)
}

wdOnsimulatedData <- function(scenario,simulationsPath){
	output <- NULL
	plotTable <- list()

	setwd(paste0(simulationsPath,scenario,'/'))
	listDaf <- lapply(list.files(pattern = 'daf',recursive=T),fread)
	listDiv <- lapply(list.files(pattern = 'div',recursive=T),fread,sep='\t',header=T)
	
	for(simulation in 1:length(listDaf)){

		# Format data from each replicate
		colnames(listDaf[[simulation]]) <- c('daf','Pi','P0')
		
		if(scenario=='length1e8'){
			listDiv[[simulation]]$m0<-1e8
			listDiv[[simulation]]$mi <- 1e8
		}
		else if(scenario=='length1e6'){
			listDiv[[simulation]]$m0<-1e6
			listDiv[[simulation]]$mi<-1e6
		}
		else{
			listDiv[[simulation]]$m0<-1e7
			listDiv[[simulation]]$mi<-1e7
		}
		
		colnames(listDiv[[simulation]]) <- c('Di', 'D0','trueAlpha','m0','mi')
		trueAlpha <- listDiv[[simulation]][,c('trueAlpha')]
		divergence <- listDiv[[simulation]][,c('Di','D0','m0','mi')]

		##DGRP
		resultDGRP0.05 <- DGRP(listDaf[[simulation]],divergence,listCutoffs = c(0.025))
		alphaDGRP0.05 <- resultDGRP0.05$Results$alpha.symbol

		resultDGRP0.1 <- DGRP(listDaf[[simulation]],divergence,listCutoffs = c(0.075))
		alphaDGRP0.1 <- resultDGRP0.1$Results$alpha.symbol
		
		alphaAsymptotic <- tryCatch({
			resultiMK <- iMKT(listDaf[[simulation]],divergence,xlow = 0, xhigh = 1)
			alphaAsymptotic <- resultiMK$`Asymptotic MK table`$alpha_asymptotic
			},error=function(e){ alphaAsymptotic0.1 <- NA; alphaAsymptotic <- NA})

		result <- data.frame(simulation,alphaDGRP0.05,alphaDGRP0.1,alphaAsymptotic,alphaAsymptotic0.1,trueAlpha)
		output <- rbind(output,result)
 
	}
	
	# Preparing data to plot
	dataPlot<-melt(output,id=c('simulation'))

	# Errors
	standardError <- mean(abs(output[['alphaStandard']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRPError <- mean(abs(output[['alphaDGRP0.05']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRP1Error <- mean(abs(output[['alphaDGRP0.1']]-output[['trueAlpha']]),na.rm=TRUE)
	FWWError <- mean(abs(output[['alphaFWW0.05']]-output[['trueAlpha']]),na.rm=TRUE)
	FWW1Error <- mean(abs(output[['alphaFWW0.1']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptoticError <- mean(abs(output[['alphaAsymptotic']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptotic1Error <- mean(abs(output[['alphaAsymptotic0.1']]-output[['trueAlpha']]),na.rm=TRUE)

	error <- matrix(c(standardError,DGRPError,DGRP1Error,FWWError,FWW1Error,asymptoticError,asymptotic1Error,0),nrow=1,ncol=8)

	meanSd <- sapply(output, function(cl) list(means=mean(cl,na.rm=TRUE), sds=sd(cl,na.rm=TRUE)))
	meanSd <- meanSd[,-1]
	meanSd <- rbind(meanSd,error)
	meanSd <- as.data.frame(meanSd)
	rownames(meanSd) <- c('mean','sd','error')

	meanSd$scenario <- scenario

	# Ploting
	plotAlpha <- ggplot(dataPlot, aes(x=variable, y=value, fill=variable)) + 
		geom_boxplot(color="grey20",alpha=0.7) + 
		labs(x = "MKT methods", y=expression(italic(α))) + 
		themePublication() + 
		scaleFillPublication(name="Method", labels=c("alphaStandard"="Standard", "alphaDGRP0.05" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW0.05" = "FWW 5%", "alphaFWW0.1" = "FWW 10%","alphaAsymptotic"="Asymptotic MKT","alphaAsymptotic0.1"="Asymptotic MKT", "trueAlpha"="Real alpha")) + 
		scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + 
		scale_x_discrete(labels=c("alphaStandard"="Standard", "alphaDGRP0.05" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW0.05" = "FWW 5%", "alphaFWW0.1" = "FWW 10%","alphaAsymptotic"="Asymptotic MKT","alphaAsymptotic0.1"="Asymptotic MKT", "trueAlpha"="Real alpha")) +  
		guides(fill=FALSE) + 
		ggtitle(paste0(scenario)) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),axis.text.y= element_text(size=20),plot.title=element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24))

	method <- rep(c("Standard","eMKT cutoff 5%","eMKT DGRP 10%","FWW cutoff 5%","MKT FWW cutoff 10%","Asymptotic MKT","Asymptotic MKT","Real alpha"),2)
					 

	plotTable[['plot']] <- plotAlpha
	plotTable[['table']] <- meanSd
	plotTable[['data']] <- dataPlot
					 
	return(plotTable)
}