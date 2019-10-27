library(splitstackshape)

samplingGenes = function(geneList,B,bins,seed=213,path=NULL){
	set.seed(seed)
	
	bootSamples = matrix(sample(geneList,size = B*bins,replace=T),B,bins)
	
	output = list()
	for (i in 1:B) {
		output[[i]] = bootSamples[i,]
	}

	if(!is.null(path)){
		save(output,file=paste0(path,'bootSampleBin',b,'.RData'))
	}

	return(output)
}

sampleAnalysis = function(data,sampling,bins,population,recomb=FALSE,xlow=0.1,xhigh=0.9){
	data = data %>% as.data.table
	data = data[pop==population]
	data$id = as.character(data$id)
	tmp = list()

	for (i in 1:length(sampling)) {
	# for (i in 2600:3500) {
		print(i)	
		sampling[[i]]           = sampling[[i]] %>% as.data.frame
		colnames(sampling[[i]]) = 'id'
		# subsetData            = data[match(data$Name, sampling[[i]]),] %>% na.omit()
		subsetData              = merge(sampling[[i]],data,by='id',all.x=T)
		
		divergence              = subsetData[,c('mi', 'di', 'm0', 'd0')]
		colnames(divergence)    = c('mi', 'Di', 'm0', 'D0')
		divergence              = divergence %>% summarize(mi=sum(mi),Di=sum(Di),m0=sum(m0),D0=sum(D0))
		
		sfs0f                   = as.data.frame(subsetData[['daf0f']])
		colnames(sfs0f)         = 'sfs0f'
		sfs0f                   = cSplit(sfs0f, 'sfs0f', ';') %>% colSums
		sfs0f                   = as.numeric(sfs0f)
		
		sfs4f                   = as.data.frame(subsetData[['daf4f']]) 
		colnames(sfs4f)         = 'sfs4f'
		sfs4f                   = cSplit(sfs4f, 'sfs4f', ';') %>% colSums
		sfs4f                   = as.numeric(sfs4f)
		
		daf                     = cbind(daf=c(0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975), sfs0f, sfs4f) %>% as.data.frame
		colnames(daf)           = c('daf','Pi','P0')
		Pi                      = sum(daf$Pi)
		P0                      = sum(daf$P0)
		D0                      = divergence$D0
		Di                      = divergence$Di
		m0                      = divergence$m0
		mi                      = divergence$mi


		if (P0 == 0 | D0  == 0 | Pi == 0 | Di == 0 | mi  == 0 | m0 == 0 | sum(m0) < (sum(P0)+sum(D0))) {
			tmp[[i]] = c(NA,NA,NA,NA,NA,NA)
			# tmp = data.frame(NA,NA,NA,NA,NA,NA)
			#output = rbind(output,tmp)
			# colnames(tmp) = c('alphaStandard','alphaEmkt1','alphaEmkt2','alphaFww1','alphaFww2','alphaAsymptotic1')
			# colnames(tmp) = c('alphaStandard','alphaEmkt1','alphaEmkt2','alphaFww1','alphaFww2','alphaAsymptotic1')
			
		}else{

			resultStandard = standardMKT(daf,divergence)
			alphaStandard = resultStandard$alpha$alpha

			PiMinus     = sum(daf[daf[['daf']] <= 0.075,'Pi'])
			PiGreater   = sum(daf[daf[['daf']] > 0.075,'Pi'])
			P0Minus     = sum(daf[daf[['daf']] <= 0.075,'P0'])
			P0Greater   = sum(daf[daf[['daf']] > 0.075,'P0'])

			if(PiGreater== 0 | P0Greater  == 0){

				alphaEmkt1 = NA
				alphaEmkt2 = NA
			}
			else{
				resultEmkt1 = eMKT(daf,divergence,listCutoffs = c(0.025))
				alphaEmkt1 = resultEmkt1$alphaCorrected$alphaCorrected
				resultEmkt2 = eMKT(daf,divergence,listCutoffs = c(0.075))
				alphaEmkt2 = resultEmkt1$alphaCorrected$alphaCorrected

			}

			resultFww1 = FWW(daf,divergence,listCutoffs = c(0.025))
			alphaFww1 = resultFww1$alphaCorrected$alphaCorrected

			resultFww2 = FWW(daf,divergence,listCutoffs = c(0.075))
			alphaFww2 = resultFww2$alphaCorrected$alphaCorrected

			##Asymptotic
			alphaAsymptotic1 = tryCatch({
				resultiMK1 = aMKT(daf,divergence,xlow = xlow, xhigh = xhigh, plot=F)
				alphaAsymptotic1 = resultiMK1$alphaCorrected$alphaAsymptotic
			},error=function(e){alphaAsymptotic1 = NA})

			# # daf1 = daf
			# # daf1$daf10 = sort(rep(seq(0.05, 0.95, 0.1), 
			# # 2))
			# # daf1 = daf1[c("daf10", "Pi", "P0")]
			# # daf1 = aggregate(. ~ daf10, data = daf1, FUN = sum)
			# # colnames(daf1) = c("daf", "Pi", "P0")
			

			tmp[[i]] = c(alphaStandard,alphaEmkt1,alphaEmkt2,alphaFww1,alphaFww2,alphaAsymptotic1)
			# tmp = data.frame(alphaStandard,alphaEmkt1,alphaEmkt2,alphaFww1,alphaFww2,alphaAsymptotic1)
			
		}
		# tmp = data.frame(alphaStandard,alphaEmkt1,alphaEmkt2,alphaFww1,alphaFww2,alphaAsymptotic,alphaAsymptotic0.1)
	}

	output = do.call(rbind.data.frame,tmp)
	output[['bin']] = bins
	colnames(output) <- c('alphaStandard','alphaEmkt1','alphaEmkt2','alphaFww1','alphaFww2','alphaAsymptotic1','bins')
	return(output)
}

wdOnConcatenation = function(sizes,pop,data){

	out = list()

	for(s in sizes){

		print(s)

		genes = unique(data[['id']]) %>% as.matrix()
		set.seed(13753)
		subsetGenes = sample(genes,size=s)

		bootSamplesBin = samplingGenes(geneList = subsetGenes,B = 3500,bins = s)
		data = data; sampling = bootSamplesBin;bins = s; population = pop

		data = data %>% as.data.table
		data = data[pop==population]
		data[['id']] = as.character(data[['id']])
		fractions = data.table()
		ci = data.table()
		xlow = 0; xhigh = 1

		for(i in 1:100){
			sampling[[i]] = sampling[[i]] %>% as.data.frame
			colnames(sampling[[i]]) = 'id'
			# subsetData = data[match(data$Name, sampling[[i]]),] %>% na.omit()
			subsetData = merge(sampling[[i]],data,by='id',all.x=T)

			divergence = subsetData[,c('mi', 'di', 'm0', 'd0')]
			colnames(divergence) = c('mi', 'Di', 'm0', 'D0')
			divergence = divergence %>% summarize(mi=sum(mi),Di=sum(Di),m0=sum(m0),D0=sum(D0))

			sfs0f = as.data.frame(subsetData[['daf0f']])
			colnames(sfs0f) = 'sfs0f'
			sfs0f = cSplit(sfs0f, 'sfs0f', ';') %>% colSums
			sfs0f = as.numeric(sfs0f)

			sfs4f = as.data.frame(subsetData[['daf4f']]) 
			colnames(sfs4f) = 'sfs4f'
			sfs4f = cSplit(sfs4f, 'sfs4f', ';') %>% colSums
			sfs4f = as.numeric(sfs4f)

			daf = cbind(daf=seq(0.05,1,0.05), sfs0f, sfs4f) %>% as.data.frame
			colnames(daf) = c('daf','Pi','P0')
			# daf=daf[daf$daf <= 0.9, ]
			D0 = divergence$D0
			Di = divergence$Di
			Pi = sum(daf$Pi)
			P0 = sum(daf$P0)

			resultsDGRP1 = DGRP(daf=daf,div=divergence,listCutoffs=0.05)
			fractions = rbind(fractions,data.table(resultsDGRP1[['Fractions']][2,],resultsDGRP1[['Fractions']][3,],'eMKT 5%'))
			resultsDGRP2 = DGRP(daf=daf,div=divergence,listCutoffs=0.1)
			fractions = rbind(fractions,data.table(resultsDGRP2[['Fractions']][2,],resultsDGRP2[['Fractions']][3,],'eMKT 10%'))
			
			
			resultsiMKT = tryCatch(iMKT(daf=daf,div=divergence,xlow=xlow,xhigh=xhigh,plot=T),error=function(e){resultiMKT=NULL})
			
			print(resultsiMKT$Fractions)
			if(is.null(resultsiMKT)){

				fractions = rbind(fractions,data.table(NA,NA,'aMKT'))
				fractions = rbind(fractions,data.table(NA,NA,'aMKT 5%'))
				fractions = rbind(fractions,data.table(NA,NA,'aMKT S'))
				fractions = rbind(fractions,data.table(NA,NA,'aMKT notCorrected'))
				
			}
			else{

				fractions = rbind(fractions,data.table(NA,resultsiMKT[['Fractions']][2,2],'aMKT'))
				fractions = rbind(fractions,data.table(NA,resultsiMKT[['Fractions']][3,2],'aMKT S'))
				fractions = rbind(fractions,data.table(NA,resultsiMKT[['Fractions']][4,2],'aMKT 5%'))
				fractions = rbind(fractions,data.table(NA,resultsiMKT[['Fractions']][5,2],'aMKT notCorrected'))
				

				ci = rbind(ci,data.table(resultsiMKT[['asymptoticMkTable']][['alpha_asymptotic']],resultsiMKT[['asymptoticMkTable']][['CI_low']],resultsiMKT[['asymptoticMkTable']][['CI_high']],as.numeric(s)))
			}
		}
		
		colnames(fractions) = c('f','b','test')
		fractions[['test']] = factor(fractions[['test']],levels=c('eMKT 5%','eMKT 10%','aMKT','aMKT S','aMKT 5%','aMKT notCorrected'))
	 	
		# df = melt(fractions)
		df = melt(fractions[,2:3])

		p = ggplot(df) + geom_boxplot(aes(x=test,y=value,fill=variable)) + themePublication() + scaleFillPublication() + xlab('Method') + ylab('Proportion of mutations') + guides(fill=guide_legend(title='Mutation type'))

		out[[paste0(s)]][['plot']] = p
		colnames(ci) = c('alpha','ciLow','ciHigh','subetLength')
		out[[paste0(s)]][['ci']] = ci
		out[[paste0(s)]][['Fractions']] = fractions

	}

	return(out)
}