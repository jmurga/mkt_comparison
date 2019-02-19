createCombinations <- function(listLengthToExpand,ratioD,ratioP){

	combinations <- list()
	
	combinationsD <- expand.grid(rep(list(1:listLengthToExpand), 2)) %>% as.data.table()
	colnames(combinationsD) <- c("Di","D0")
	combinationsD$DivRatio <- combinationsD$D0/combinationsD$Di
	combinationsD <- combinationsD[DivRatio==ratioD]

	combinationsD$ND <- combinationsD$Di+combinationsD$D0


	combinationsP<-expand.grid(rep(list(1:listLengthToExpand), 2)) %>% as.data.table()
	colnames(combinationsP) <- c("Pi","P0")
	combinationsP$PolRatio <- combinationsP$P0/combinationsP$Pi

	combinationsP <- combinationsP[PolRatio == ratioP]
	combinationsP$NP <- combinationsP$Pi+combinationsP$P0
	
	combinations <- cbind(combinationsD, combinationsP)
	
	
	return(combinations)
}


tableToIter <- function(data=NULL,iterations,combinations=NULL) {
	if(combinations == 'divergence'){
		output <- list()
		output[['results']] <- data.frame()
		plots <- list()

		for (div in 1:length(data$Di)) {
			# print(div)
			di <- data[div,]$Di
			d0 <- data[div,]$D0

			tmp <- NULL

			for (i in 1:iterations) {
				pi <- i
				p0 <- data[div,]$D0

				mktTable <- data.frame(Polymorphism = c(p0,pi), Divergence = c(d0, di), row.names = c("Neutral class", "Selected class"))
				alpha <- 1 - ((pi*d0)/(p0*di))				
				pvalue <- fisher.test(mktTable)$p.value
				
				results <- data.frame(pi,p0,di,d0,alpha,pvalue)
				tmp <- rbind(tmp,results)
		
			}

			tmpSubset <- cbind(di,d0, ">Pi"=tail(subset(tmp,tmp$pvalue<0.05 & tmp$alpha>0),1)[1,1],sign_pos=sum(tmp$pvalue<0.05 & tmp$alpha>0),">Pi"=head(subset(tmp,tmp$pvalue<0.05 & tmp$alpha<0),1)[1,1], sign_neg=sum(tmp$pvalue<0.05 & tmp$alpha<0))

				p <- ggplot(tmp, aes(x = pi)) + geom_point(aes(y = alpha,  colour = pvalue < 0.05)) + 
				scale_y_continuous(breaks = pretty(tmp$alpha, n = 5)) +
				scaleColourPublication(name="P-value") +
				labs(y = expression(italic(α)),
					 x = "Pi") + 
				theme(legend.position = c(0.8, 0.9))  + themePublication() + ggtitle(paste("Di: ", di, " - D0: ", d0))
			  	
			  
			  	plots[[div]] <- p
				output[['results']] <- rbind(output[['results']],tmpSubset)
		}
		output[['plots']] <- plot_grid(plotlist=plots,ncol=6)
	}
	else if(combinations == 'polymorphism'){

		output <- list()
		plots <- list()
		for (div in 1:length(data$Di)) {
			# print(div)

			di <- data[div,]$Di
			d0 <- data[div,]$D0
			nd <- data[div,]$ND
			
			tmp <- NULL
			
			for (i in 1:length(data$Pi)) {
			 
				pi <- data[i,]$Pi
				p0 <- data[i,]$P0
				np <- data[i,]$NP
				
				mktTable <- data.frame(P = c(p0, pi),Divergence = c(d0, di), row.names = c("Neutral class", "Selected class"))
				alpha <- 1 - ((di*p0)/(p0*di))
				pvalue <- fisher.test(mktTable)$p.value
				
				result <- data.frame(pi,p0,di,d0,alpha,pvalue,nd,np)
				tmp<-rbind(tmp,result)
			}

			output[['results']] <- rbind(output[['results']],tmp)
			ratioP <- output[['results']]$p0[1]/output[['results']]$pi[1]
			ratioD <- output[['results']]$di[1]/output[['results']]$p0[1]

			p<-ggplot(output[['results']])
			p <- p+geom_tile(aes(np,nd,fill=pvalue))+ggtitle(paste0("α = ",unique(alpha))) + xlab(paste0("Total polymorphic sites (",ratioP,"P0/Pi)")) + ylab(paste0("Total diverge sites (D0/",ratioD,"Di)"))  + scale_fill_distiller(palette = "Spectral",name="P-value")+theme(axis.text.x = element_text(size=14, angle=90)) + themePublication() + scale_x_continuous(expand = c(0,0),breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(expand = c(0,0),breaks = scales::pretty_breaks(n = 20))
			
			output[['plot']] <- p
			# data <- do.call(rbind.data.frame, output)


			}
	}
	else if(is.null(combinations)){

		di <- round(mean(data$di))
		d0 <- round(mean(data$d0))

		pi<-di
		p0<-d0

		output <- list()
		tmp <- NULL

		for (i in seq(0,100,1)) {
		  P0<-D0
		  Pi<-i
		  mkt_table <- data.frame(Polymorphism = c(p0, pi), Divergence = c(d0, di), row.names = c("Neutral class","Selected class"))
		  alpha <- 1 - (mkt_table[2, 1]/mkt_table[1, 1]) * (mkt_table[1,2]/mkt_table[2, 2])
		  pvalue <- fisher.test(mkt_table)$p.value
		  
		  result <- data.frame(i,alpha,pvalue)
		  tmp <- rbind(tmp,result)

		}

		p <- ggplot(tmp, aes(x = i)) + geom_point(aes(y = alpha,  colour = pvalue < 0.05)) + scale_y_continuous(breaks = pretty(tmp$alpha, n = 5)) +scaleColourPublication() +labs(y = "alpha",x = "Pi") + theme(legend.position = c(0.8, 0.9))  + themePublication() 
		
		output[['results']] <- tmp
		output[['plot']] <- p
	}
	return(output)
}
