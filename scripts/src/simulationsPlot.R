library(iMKT)
library(forcats)
library(cowplot)
library(reshape2)
library(data.table)
library(dplyr)
source('/home/jmurga/mkt/201902/scripts/src/plotStyle.R')

mktOnsimulatedData <- function(scenario,simulationsPath='/home/jmurga/mkt/201902/rawData/simulations/alphaComparisons/'){
	
	output <- NULL
	plotTable <- list()

	setwd(paste0(simulationsPath,'/',scenario,'/'))
	# setwd(paste0('/home/jmurga/mkt/201902/rawData/simulations/alphaComparisons/l1e6'))
	listDaf <- lapply(list.files(pattern = 'daf',recursive=T),read.table,sep='\t',header=T)
	listDiv <- lapply(list.files(pattern = 'div',recursive=T),read.table,sep='\t',header=T)
	
	for(simulation in 1:length(listDaf)){
	# for(simulation in 1:10){

		print(simulation)

		daf <- listDaf[[simulation]][,c('daf','Pi','P0')]
		divergence <- listDiv[[simulation]][,c('Di','D0','m0','mi')]
		trueAlpha <- listDiv[[simulation]][,c('trueAlpha')]

		# Execute each mkt test through iMKT package
		## StandardMKT
		resultStandard <- standardMKT(daf,divergence)
		alphaStandard <- resultStandard$alpha$alpha

		##DGRP
		resultDGRP <- eMKT(daf,divergence, listCutoffs=c(0.025,0.075), plot=F)
		alphaDGRP <- resultDGRP$alphaCorrected$alphaCorrected
		##FWW
		resultFWW <- FWW(daf,divergence,listCutoffs = c(0.025,0.075), plot=F)
		alphaFWW <- resultFWW$alphaCorrected$alpha

		##Asymptotic
		alphaAsymptotic1 <- tryCatch({
			resultiMK1 <- aMKT(daf,divergence,xlow = 0, xhigh = 1, plot=F)
			alphaAsymptotic1 <- resultiMK1$alphaCorrected$alphaAsymptotic
			},error=function(e){alphaAsymptotic1 <- NA})

		alphaAsymptotic2 <- tryCatch({
			resultiMK2 <- aMKT(daf,divergence,xlow = 0.1, xhigh = 0.9, plot=F)
			alphaAsymptotic2 <- resultiMK2$alphaCorrected$alphaAsymptotic
			},error=function(alphaAsymptotic2){ alphaAsymptotic2 <- NA}
		)

		# result <- data.table(simulation,alphaStandard,alphaDGRP,alphaFWW1,alphaFWW2,alphaAsymptotic1,alphaAsymptotic2,trueAlpha)
		result <- data.table(simulation,alphaStandard,alphaDGRP[1],alphaDGRP[2],alphaFWW[1],alphaFWW[2],alphaAsymptotic1,alphaAsymptotic2,trueAlpha)
		colnames(result) <- c('simulation','alphaStandard','alphaDGRP1','alphaDGRP2','alphaFWW1','alphaFWW2','alphaAsymptotic1','alphaAsymptotic2','trueAlpha')
		output <- rbind(output,result)
	}
	
	# Preparing data to plot
	dataPlot <- reshape2::melt(output,id=c('simulation'))

	# Errors
	standardError <- mean(abs(output[['alphaStandard']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRPError <- mean(abs(output[['alphaDGRP1']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRP2Error <- mean(abs(output[['alphaDGRP2']]-output[['trueAlpha']]),na.rm=TRUE)
	FWWError <- mean(abs(output[['alphaFWW1']]-output[['trueAlpha']]),na.rm=TRUE)
	FWW1Error <- mean(abs(output[['alphaFWW2']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptoticError <- mean(abs(output[['alphaAsymptotic1']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptotic1Error <- mean(abs(output[['alphaAsymptotic2']]-output[['trueAlpha']]),na.rm=TRUE)

	error <- matrix(c(standardError,DGRPError,DGRP2Error,FWWError,FWW1Error,asymptoticError,asymptotic1Error,0),nrow=1,ncol=8)

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
		scaleFillPublication(name="Method", labels=c("alphaStandard"="Standard", "alphaDGRP1" = "eMKT 5%", "alphaDGRP2" = "eMKT 10%", "alphaFWW1" = "FWW 5%", "alphaFWW2" = "FWW 10%","alphaAsymptotic1"="Asymptotic MKT","alphaAsymptotic2"="Asymptotic MKT", "trueAlpha"="True alpha")) + 
		scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + 
		scale_x_discrete(labels=c("alphaStandard"="Standard", "alphaDGRP1" = "eMKT 5%", "alphaDGRP2" = "eMKT 10%", "alphaFWW1" = "FWW 5%", "alphaFWW2" = "FWW 10%","alphaAsymptotic1"="Asymptotic MKT","alphaAsymptotic2"="Asymptotic MKT 0.1-0.9", "trueAlpha"="True alpha")) +  
		guides(fill=FALSE) + 
		ggtitle(paste0(scenario)) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),axis.text.y= element_text(size=20),plot.title=element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24))

	# method <- rep(c("Standard","eMKT","FWW cutoff 5%","MKT FWW cutoff 10%","Asymptotic MKT","Asymptotic MKT 0.1-0.9","True alpha"),2)
					 

	plotTable[['plot']] <- plotAlpha
	plotTable[['table']] <- meanSd
	plotTable[['data']] <- dataPlot
					 
	return(plotTable)
}

wdOnsimulatedData <- function(scenario,simulationsPath){

	output <- NULL
	plotTable <- list()

	setwd(paste0(simulationsPath,'/',scenario,'/'))
	# setwd(paste0('/home/jmurga/mkt/201902/rawData/simulations/alphaComparisons/l1e6'))
	listDaf <- lapply(list.files(pattern = 'daf',recursive=T),read.table,sep='\t',header=T)
	listDiv <- lapply(list.files(pattern = 'div',recursive=T),read.table,sep='\t',header=T)
	
	for(simulation in 1:length(listDaf)){
	# for(simulation in 1:10){

		print(simulation)

		daf <- listDaf[[simulation]][,c('daf','Pi','P0')]
		divergence <- listDiv[[simulation]][,c('Di','D0','m0','mi')]
		trueAlpha <- listDiv[[simulation]][,c('trueAlpha')]
		trueB <- listDiv[[simulation]][,c('b')]

		##DGRP
		resultDGRP <- eMKT(daf,divergence, listCutoffs=c(0.025,0.075), plot=F)
		bDGRP <- resultDGRP$fractions$b

		##Asymptotic
		bAsymptotic1 <- tryCatch({
			resultiMK1 <- aMKT(daf,divergence,xlow = 0, xhigh = 1, plot=F)
			bAsymptotic1 <- resultiMK1$fractions$b
			},
			error=function(e)
			{
				bAsymptotic1 <- NA
			}
		)

		# if((resultiMK1$alphaCorrected$ciHigh - resultiMK1$alphaCorrected$ciLow) > 1){
		# 	bAsymptotic1 <- NA
		# }

		bAsymptotic2 <- tryCatch({
			resultiMK2 <- aMKT(daf,divergence,xlow = 0.1, xhigh = 0.9, plot=F)
			bAsymptotic2<- resultiMK2$fractions$b
			},
			error=function(e){ 
				bAsymptotic2 <- NA
			}
		)

		# if((resultiMK2$alphaCorrected$ciHigh - resultiMK2$alphaCorrected$ciLow) > 1){
		# 	bAsymptotic1 <- NA
		# }

		result <- data.table(simulation,bDGRP[1],bDGRP[2],bAsymptotic1,bAsymptotic2,trueB)
		colnames(result) <- c('simulation','bDGRP1','bDGRP2','bAsymptotic1','bAsymptotic2','trueB')
		output <- rbind(output,result)
	}
	
	# Preparing data to plot
	dataPlot <- reshape2::melt(output,id=c('simulation'))

	# Errors
	DGRPError <- mean(abs(output[['bDGRP1']]-output[['trueB']]),na.rm=TRUE)
	DGRP2Error <- mean(abs(output[['bDGRP2']]-output[['trueB']]),na.rm=TRUE)
	asymptoticError <- mean(abs(output[['bAsymptotic1']]-output[['trueB']]),na.rm=TRUE)
	asymptotic1Error <- mean(abs(output[['bAsymptotic2']]-output[['trueB']]),na.rm=TRUE)

	error <- matrix(c(DGRPError,DGRP2Error,asymptoticError,asymptotic1Error,0),nrow=1,ncol=5)

	meanSd <- sapply(output, function(cl) list(means=mean(cl,na.rm=TRUE), sds=sd(cl,na.rm=TRUE)))
	meanSd <- meanSd[,-1]
	meanSd <- rbind(meanSd,error)
	meanSd <- as.data.frame(meanSd)
	rownames(meanSd) <- c('mean','sd','error')

	meanSd$scenario <- scenario

	# Ploting
	plotAlpha <- ggplot(dataPlot, aes(x=variable, y=value, fill=variable)) + 
		geom_boxplot(color="grey20",alpha=0.7) + 
		labs(x = "MKT methods", y='b') + 
		themePublication() + 
		scaleFillPublication(name="Method", labels=c('bDGRP1' = 'eMKT 5%', 'bDGRP2' = 'eMKT 10%', 'bFWW1' = 'FWW 5%', 'bFWW2' = 'FWW 10%','bAsymptotic1'='Asymptotic MKT','bAsymptotic2'='Asymptotic MKT', 'trueB'='True b')) + 
		scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + 
		scale_x_discrete(labels=c('bStandard'='Standard', 'bDGRP1' = 'eMKT 5%', 'bDGRP2' = 'eMKT 10%','bAsymptotic1'='Asymptotic MKT','bAsymptotic2'='Asymptotic MKT 0.1-0.9', 'trueB'='True b')) +  
		guides(fill=FALSE) + 
		ggtitle(paste0(scenario)) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),axis.text.y= element_text(size=20),plot.title=element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24))

	# method <- rep(c("Standard","eMKT","FWW cutoff 5%","MKT FWW cutoff 10%","Asymptotic MKT","Asymptotic MKT 0.1-0.9","True alpha"),2)
					 

	plotTable[['plot']] <- plotAlpha
	plotTable[['table']] <- meanSd
	plotTable[['data']] <- dataPlot
					 
	return(plotTable)
					 
}