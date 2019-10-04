library(iMKT)
library(forcats)
library(cowplot)
library(reshape2)
library(data.table)
source('/home/jmurga/mkt/201902/scripts/src/plotStyle.R')

mktOnsimulatedData <- function(scenario,simulationsPath){
	
	output <- NULL
	plotTable <- list()

	setwd(paste0(simulationsPath,'/',scenario,'/'))
	# setwd(paste0('/home/jmurga/mkt/201902/rawData/simulations/alphaComparisons/l1e6'))
	listDaf <- lapply(list.files(pattern = 'daf',recursive=T),fread)
	listDiv <- lapply(list.files(pattern = 'div',recursive=T),fread,sep='\t',header=T)
	
	for(simulation in 1:length(listDaf)){

		print(simulation)

		daf <- listDaf[[simulation]][,c('daf','Pi','P0')]
		divergence <- listDiv[[simulation]][,c('Di','D0','m0','mi')]
		trueAlpha <- listDiv[[simulation]][,c('trueAlpha')]

		# Execute each mkt test through iMKT package
		## StandardMKT
		resultStandard <- standardMKT(daf,divergence)
		alphaStandard <- resultStandard$alpha$alpha

		##DGRP
		resultDGRP <- eMKT(daf,divergence, plot=F)
		alphaDGRP <- resultDGRP$alphaCorrected$alphaCorrected

		##FWW
		resultFWW1 <- FWW(daf,divergence,listCutoffs = c(0.025), plot=F)
		alphaFWW1 <- resultFWW1$alphaCorrected$alpha

		resultFWW2 <- FWW(daf,divergence,listCutoffs = c(0.075), plot=F)
		alphaFWW2 <- resultFWW2$alphaCorrected$alpha

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

		result <- data.table(simulation,alphaStandard,alphaDGRP,alphaFWW1,alphaFWW2,alphaAsymptotic1,alphaAsymptotic2,trueAlpha)
		output <- rbind(output,result)
	}
	
	# Preparing data to plot
	dataPlot <- reshape2::melt(output,id=c('simulation'))

	# Errors
	standardError <- mean(abs(output[['alphaStandard']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRPError <- mean(abs(output[['alphaDGRP']]-output[['trueAlpha']]),na.rm=TRUE)
	FWWError <- mean(abs(output[['alphaFWW1']]-output[['trueAlpha']]),na.rm=TRUE)
	FWW1Error <- mean(abs(output[['alphaFWW2']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptoticError <- mean(abs(output[['alphaAsymptotic1']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptotic1Error <- mean(abs(output[['alphaAsymptotic2']]-output[['trueAlpha']]),na.rm=TRUE)

	error <- matrix(c(standardError,DGRPError,FWWError,FWW1Error,asymptoticError,asymptotic1Error,0),nrow=1,ncol=7)

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
		scaleFillPublication(name="Method", labels=c("alphaStandard"="Standard", "alphaDGRP" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW1" = "FWW 5%", "alphaFWW2" = "FWW 10%","alphaAsymptotic1"="Asymptotic MKT","alphaAsymptotic2"="Asymptotic MKT", "trueAlpha"="True alpha")) + 
		scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + 
		scale_x_discrete(labels=c("alphaStandard"="Standard", "alphaDGRP" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW1" = "FWW 5%", "alphaFWW2" = "FWW 10%","alphaAsymptotic1"="Asymptotic MKT","alphaAsymptotic2"="Asymptotic MKT 0.1-0.9", "trueAlpha"="True alpha")) +  
		guides(fill=FALSE) + 
		ggtitle(paste0(scenario)) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),axis.text.y= element_text(size=20),plot.title=element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24))

	method <- rep(c("Standard","eMKT","FWW cutoff 5%","MKT FWW cutoff 10%","Asymptotic MKT","Asymptotic MKT 0.1-0.9","True alpha"),2)
					 

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
		colnames(daf) <- c('daf','Pi','P0')
		
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
		resultDGRP <- DGRP(daf,divergence,listCutoffs = c(0.025))
		alphaDGRP <- resultDGRP0.05$Results$alpha.symbol

		resultDGRP0.1 <- DGRP(daf,divergence,listCutoffs = c(0.075))
		alphaDGRP0.1 <- resultDGRP0.1$Results$alpha.symbol
		
		alphaAsymptotic1 <- tryCatch({
			resultiMK1 <- aMKT(daf,divergence,xlow = 0, xhigh = 1)
			alphaAsymptotic1 <- resultiMK1$`Asymptotic MK table`$alpha_asymptotic
			},error=function(e){ alphaAsymptotic2 <- NA; alphaAsymptotic1 <- NA})

		result <- data.frame(simulation,alphaDGRP,alphaDGRP0.1,alphaAsymptotic1,alphaAsymptotic2,trueAlpha)
		output <- rbind(output,result)
 
	}
	
	# Preparing data to plot
	dataPlot<-melt(output,id=c('simulation'))

	# Errors
	standardError <- mean(abs(output[['alphaStandard']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRPError <- mean(abs(output[['alphaDGRP']]-output[['trueAlpha']]),na.rm=TRUE)
	DGRP1Error <- mean(abs(output[['alphaDGRP0.1']]-output[['trueAlpha']]),na.rm=TRUE)
	FWWError <- mean(abs(output[['alphaFWW1']]-output[['trueAlpha']]),na.rm=TRUE)
	FWW1Error <- mean(abs(output[['alphaFWW2']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptoticError <- mean(abs(output[['alphaAsymptotic1']]-output[['trueAlpha']]),na.rm=TRUE)
	asymptotic1Error <- mean(abs(output[['alphaAsymptotic2']]-output[['trueAlpha']]),na.rm=TRUE)

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
		scaleFillPublication(name="Method", labels=c("alphaStandard"="Standard", "alphaDGRP" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW1" = "FWW 5%", "alphaFWW2" = "FWW 10%","alphaAsymptotic1"="Asymptotic MKT","alphaAsymptotic2"="Asymptotic MKT", "trueAlpha"="Real alpha")) + 
		scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + 
		scale_x_discrete(labels=c("alphaStandard"="Standard", "alphaDGRP" = "eMKT 5%", "alphaDGRP0.1" = "eMKT 10%", "alphaFWW1" = "FWW 5%", "alphaFWW2" = "FWW 10%","alphaAsymptotic1"="Asymptotic MKT","alphaAsymptotic2"="Asymptotic MKT", "trueAlpha"="Real alpha")) +  
		guides(fill=FALSE) + 
		ggtitle(paste0(scenario)) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),axis.text.y= element_text(size=20),plot.title=element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24))

	method <- rep(c("Standard","eMKT cutoff 5%","eMKT DGRP 10%","FWW cutoff 5%","MKT FWW cutoff 10%","Asymptotic MKT","Asymptotic MKT","Real alpha"),2)
					 

	plotTable[['plot']] <- plotAlpha
	plotTable[['table']] <- meanSd
	plotTable[['data']] <- dataPlot
					 
	return(plotTable)
}