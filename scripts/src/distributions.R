distributionSampled <- function(data,sampling){
	for (pop in unique(data$Population)){

		print(pop)

		subsetData <- data[data$Population == pop, ]
		genes <- subsetData$symbol
		standardList <- rep(0,100)
		FWWList <- rep(0,100)
		DGRPList <- rep(0,100)


		for(i in 1:100){


			subsetGenes <- as.character(sample(genes,sampling,replace=T))
			standard <- NULL
			
			while(is.null(standard)){
				subsetGenes <- as.character(sample(genes,sampling,replace=F))
				try(standard <- PopHumanAnalysis(genes=subsetGenes,pops=pop,recomb=F,test='standardMKT'))
			}

			standardList[i] <- standard[[paste0('Population =  ',pop)]]$alpha.symbol

			fww <- PopHumanAnalysis(genes=subsetGenes,pops=pop,recomb=F,test='FWW')
			FWWList[i] <- fww[[paste0('Population =  ',pop)]]$Results$alpha.symbol[2]

			dgrp <- PopHumanAnalysis(genes=subsetGenes,pops=pop,recomb=F,test='DGRP')
			DGRPList[i] <- dgrp[[paste0('Population =  ',pop)]]$Results$alpha.symbol[2]
			print(DGRPList)
		}

		if(length(subsetGenes) <25){
			print('here')
			list10Standard <- standardList
			list10FWW <- FWWList
			list10DGRP <- DGRPList

			save(list10Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'10.RData'))
			save(list10FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'10.RData'))
			save(list10DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'10.RData'))

		}
		else if(length(subsetGenes) >=25 && length(subsetGenes) < 50){
			list25Standard <- standardList
			list25FWW <- FWWList
			list25DGRP <- DGRPList
			save(list25Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'25.RData'))
			save(list25FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'25.RData'))
			save(list25DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'25.RData'))
		}
		else if(length(subsetGenes) >=50 && length(subsetGenes) < 100){
			list50Standard <- standardList
			list50FWW <- FWWList
			list50DGRP <- DGRPList
			save(list50Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'50.RData'))
			save(list50FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'50.RData'))
			save(list50DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'50.RData'))
		}
		else if(length(subsetGenes) >=100 && length(subsetGenes) < 200){
			list100Standard <- standardList
			list100FWW <- FWWList
			list100DGRP <- DGRPList
			save(list100Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'100.RData'))
			save(list100FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'100.RData'))
			save(list100DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'100.RData'))
		}
		else if(length(subsetGenes) >=200 && length(subsetGenes) < 300){
			list200Standard <- standardList
			list200FWW <- FWWList
			list200DGRP <- DGRPList
			save(list200Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'200.RData'))
			save(list200FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'200.RData'))
			save(list200DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'200.RData'))
		}
		else if(length(subsetGenes) >=300 && length(subsetGenes) < 400){
			list300Standard <- standardList
			list300FWW <- FWWList
			list300DGRP <- DGRPList
			save(list300Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'300.RData'))
			save(list300FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'300.RData'))
			save(list300DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'300.RData'))
		}
		else if(length(subsetGenes) >=400 && length(subsetGenes) < 500){
			list400Standard <- standardList
			list400FWW <- FWWList
			list400DGRP <- DGRPList
			save(list400Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'400.RData'))
			save(list400FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'400.RData'))
			save(list400DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'400.RData'))
		}
		else{
			list500Standard <- standardList
			list500FWW <- FWWList
			list500DGRP <- DGRPList
			save(list500Standard,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/standard',pop,'500.RData'))
			save(list500FWW,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/fww',pop,'500.RData'))
			save(list500DGRP,file=paste0('/home/jmurga/positiveSelectionHuman/201901/rawData/humans/distributions/dgrp',pop,'500.RData'))
		}
	}
}
