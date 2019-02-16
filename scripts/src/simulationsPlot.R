library(iMKT)
library(forcats)
library(cowplot)
library(reshape2)

mktOnsimulatedData <- function(scenario,simulationsPath){
	output <-NULL
	imkOutput<-NULL
	
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

		resultDGRP0.15 <- DGRP(listDaf[[simulation]],divergence,listCutoffs = c(0.075))
		alphaDGRP0.15 <- resultDGRP0.15$Results$alpha.symbol
		
		##FWW
		resultFWW0.05 <- FWW(listDaf[[simulation]],divergence,listCutoffs = c(0.025))
		alphaFWW0.05 <- resultFWW0.05$Results$alpha.symbol

		resultFWW0.15 <- FWW(listDaf[[simulation]],divergence,listCutoffs = c(0.075))
		alphaFWW0.15 <- resultFWW0.15$Results$alpha.symbol

		##Asymptotic
		tryCatch({
			resultiMK0.1 <- iMKT(listDaf[[simulation]],divergence,xlow = 0.1, xhigh = 0.9)
			alphaAsymptotic0.1 <- resultiMK0.1$`Asymptotic MK table`$alpha_asymptotic}
			resultiMK <- iMKT(listDaf[[simulation]],divergence,xlow = 0, xhigh = 1)
			alphaAsymptotic <- resultiMK$`Asymptotic MK table`$alpha_asymptotic},error=function(e){ alphaAsymptotic0.1 <- NA;alphaAsymptotic <- NA}
		)

		result <- data.frame(simulation,alphaStandard,alphaDGRP0.05,alphaDGRP0.15,alphaFWW0.05,alphaFWW0.15,alphaAsymptotic0.1,alphaAsymptotic,trueAlpha)
		output <- rbind(output,result)
 
	}
	
	# Preparing data to plot
	dataPlot<-melt(output,id=c('simulation'))
	
	
	mean_sd <-sapply(output, function(cl) list(means=mean(cl,na.rm=TRUE), sds=sd(cl,na.rm=TRUE)))

	# Ploting
	plotAlpha <- ggplot(dataPlot, aes(x=variable, y=value, fill=variable)) + geom_boxplot(color="grey20",alpha=0.7) + labs(x = "MKT methods", y=expression(italic(α))) + themePublication() + scaleFillPublication(name="Method", labels=c("alpha_standard"="MKT standard", "alpha_DGRP0.05" = "MKT DGRP cutoff 0.05", "alpha_DGRP0.15" = "MKT DGRP cutoff 0.15", "alpha_FWW0.05" = "MKT FWW cutoff 0.05", "alpha_FWW0.15" = "MKT FWW cutoff 0.15","alpha_iMK"="iMKT", "true_alpha"="Real alpha")) + scale_y_continuous(breaks = pretty(dataPlot$value, n = 5)) + scale_x_discrete(labels=c("alpha_standard"="Standard", "alpha_DGRP0.05" = "DGRP 5%", "alpha_DGRP0.15" = "DGRP 10%", "alpha_FWW0.05" = "FWW 5%", "alpha_FWW0.15" = "FWW 10%","alpha_iMK"="iMKT", "true_alpha"="Real alpha")) +  guides(fill=FALSE) + ggtitle(paste0(scenario))

	method <- rep(c("Standard","MKT DGRP cutoff 12.5%","MKT DGRP cutoff 17.5%","FWW DGRP cutoff 12.5%","MKT FWW cutoff 17.5%","iMKT","Real alpha"),2)
	type <- sort(rep(c(T,F),7), decreasing = T)
					 
	fractions <- c(1,1,1,1,1,0.9,1,0,0,0,0,0,0.1,0)
	fraction <- data.frame(method,fractions,type)

	fractionPlot<-ggplot(fraction) + geom_bar(stat="identity", aes_string(x=fct_inorder(as.factor(method)), y=fractions, fill=type), color="black") + themePublication() + ylab(label="Fraction") + xlab(label="Bins")  +     scaleFillPublication(name="Number of analyzable simulations") + theme(axis.line=element_blank())
					 
	plot <- plot_grid(plotAlpha + theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()),fractionPlot,ncol = 1,rel_heights = c(1,0.25))

	dataPlot$type<-scenario
	
	#  Printeado por pantalla, no sabemos si se guardó
	# table(dataPlot$variable)
	# mean_sd
	# quantile(output_baseline$alpha_standard,c(0.025,0.975))
	# quantile(output_baseline$alpha_DGRP0.05,c(0.025,0.975))   
	# quantile(output_baseline$alpha_DGRP0.15,c(0.025,0.975))   
	# quantile(output_baseline$alpha_FWW0.05,c(0.025,0.975))   
	# quantile(output_baseline$alpha_FWW0.15,c(0.025,0.975))   
	# quantile(output_baseline$alpha_iMK,c(0.025,0.975))
	# quantile(output_baseline$true_alpha,c(0.025,0.975))
					 
	return(plotAlpha)

}