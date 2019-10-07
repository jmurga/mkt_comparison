library(iMKT)
library(forcats)
library(cowplot)
library(reshape2)
library(data.table)
library(dplyr)
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

		print(simulation)

		daf <- listDaf[[simulation]][,c('daf','Pi','P0')]
		divergence <- listDiv[[simulation]][,c('Di','D0','m0','mi')]
		trueAlpha <- listDiv[[simulation]]$trueAlpha
		trueB <- listDiv[[simulation]]$b

		# Execute each mkt test through iMKT package
		##DGRP
		resultDGRP <- eMKT(daf,divergence, plot=F)
		alphaDGRP <- resultDGRP$alphaCorrected$alphaCorrected
		bDGRP <- resultDGRP$fractions$b

		##Asymptotic
		resultiMK1 <- tryCatch({
			aMKT(daf,divergence,xlow = 0, xhigh = 1, plot=F)
			},error=function(e){
				return(NULL)
			}
		)

		resultiMK2 <- tryCatch({
			aMKT(daf,divergence,xlow = 0.1, xhigh = 0.9, plot=F)
			},error=function(e){ 
				return(NULL)
			}
		)

		if(is.null(resultiMK1))
		{
			alphaAsymptotic1 <- NA
			alphaAsymptotic2 <- NA
			bAsymptotic1 <- NA
			bAsymptotic2 <- NA
		}
		else if(is.null(resultiMK2))
		{
			alphaAsymptotic1 <- resultiMK1$alphaCorrected$alphaAsymptotic
			bAsymptotic1 <- resultiMK1$Fractions[3,1]
			alphaAsymptotic2 <- NA
			bAsymptotic2 <- NA
		}
		else
		{
			alphaAsymptotic1 <- resultiMK1$alphaCorrected$alphaAsymptotic
			bAsymptotic1 <- resultiMK1$Fractions[3,1]
			alphaAsymptotic2 <- resultiMK2$alphaCorrected$alphaAsymptotic
			bAsymptotic2 <- resultiMK2$Fractions[3,1]
		}

		resultAlpha <- data.table('DGRP'=alphaDGRP,'asymp1'=alphaAsymptotic1,'asymp2'=alphaAsymptotic2,'trueValue'=trueAlpha,'measure'='alpha')
		resultB <- data.table('DGRP'=bDGRP,'asymp1'=bAsymptotic1,'asymp2'=bAsymptotic2,'trueValue'=trueB,'measure'='b')

		result <- rbind(resultAlpha,resultB)
		output <- rbind(output,result)
	}
	
	# Preparing data to plot
	dataPlot <- reshape2::melt(output,id=c('measure'))

	# Errors
	errorB <- output[output$measure == 'b',]
	DGRPError <- mean(abs(errorB[['DGRP']]-errorB[['trueValue']]),na.rm=TRUE)
	asymptotic1Error <- mean(abs(errorB[['asymp1']]-errorB[['trueValue']]),na.rm=TRUE)
	asymptotic2Error <- mean(abs(errorB[['asymp2']]-errorB[['trueValue']]),na.rm=TRUE)

	error <- matrix(c(DGRPError,asymptotic1Error,asymptotic2Error,0),nrow=1,ncol=4)

	meanSd <- sapply(errorB[,1:4], function(cl) list(means=mean(cl,na.rm=TRUE), sds=sd(cl,na.rm=TRUE)))
	meanSd <- rbind(meanSd,error)
	meanSd <- as.data.frame(meanSd)
	rownames(meanSd) <- c('mean','sd','error')

	meanSd$scenario <- scenario

	# Ploting
	dp <- dataPlot[dataPlot$measure=='b',]
	p <- ggplot(dp, aes(x=variable, y=value, fill=measure)) + 
		geom_boxplot(color="grey20",alpha=0.7,show.legend=F) + 
		labs(x = "MKT methods",y="B estimations") + 
		themePublication() + 
		scaleFillPublication(name="Method", labels=c("DGRP" = "eMKT","asymp1"="Asymptotic MKT","asymp2"="Asymptotic MKT 0.1-0.9", "trueValue"="True alpha/b")) + 
		scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + 
		scale_x_discrete(labels=c("DGRP" = "eMKT","asymp1"="Asymptotic MKT","asymp2"="Asymptotic MKT 0.1-0.9", "trueValue"="b")) +  
		ggtitle(paste0(scenario)) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),axis.text.y= element_text(size=20),plot.title=element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24))

	method <- rep(c("Standard","eMKT","FWW cutoff 5%","MKT FWW cutoff 10%","Asymptotic MKT","Asymptotic MKT 0.1-0.9","True alpha"),2)
					 

	plotTable[['plot']] <- p
	plotTable[['table']] <- meanSd
	plotTable[['data']] <- dataPlot
					 
	return(plotTable)
}