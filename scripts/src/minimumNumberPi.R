createCombinations <- function(listLengthToExpand,ratioD,ratioP){

    combinations <- list()
    
    combinationsD <- expand.grid(rep(list(1:listLengthToExpand), 2)) %>% as.data.table()
    colnames(combinationsD) <- c("Di","D0")
    combinationsD$DivRatio <- combinationsD$D0/combinationsD$Di
    combinationsD <- combinationsD[DivRatio==ratioD]

    combinationsD$ND <- combinationsD$Di+combinationsD$D0


    combinationsP<-expand.grid(rep(list(1:listLengthToExpand), 2)) %>% as.data.table()
    colnames(combinationsP) <- c("Pi","P0")
    combinationsP$DivRatio <- combinationsP$P0/combinationsP$Pi

    combinationsP <- combinationsP[DivRatio == ratioP]
    combinationsP$NP <- combinationsP$Pi+combinationsP$P0
    
    combinations <- cbind(combinationsD, combinationsP)
    combinations <- 
    
    return(combinations)
}


tableToIter <- function(data=NULL,di,d0,iterations,subset=NULL,ratio=NULL) {
    
    output <- list()
    tmp <- NULL
    
    for (i in 1:iterations) {
        tmp[i, 1]
        P0 <- d0
        Pi <- i
        if(!is.null(data)){
            Pi <- data[i,]$Pi
            P0 <- data[i,]$P0
            if(!is.null(ratio)){
                NP <- data[i,]$NP
            }
        }
        mktTable <- data.frame(Polymorphism = c(P0, Pi), Divergence = c(D0, Di), 
            row.names = c("Neutral class", "Selected class"))
        
        alpha <- 1 - (mktTable[2, 1]/mktTable[1, 1]) * (mktTable[1, 2]/mktTable[2, 
            2])
        
        pvalue <- fisher.test(mktTable)$p.value
        
        if(!is.null(ratio)){
            result <- data.frame(Pi,P0,Di,D0,alpha,pvalue,NP) 
        }
        else{
            result <- data.frame(Pi,P0,Di,D0,alpha,pvalue) 
        }
        
        tmp <- rbind(tmp, result)
    }

    if(!is.null(subset)){

        tmpSubset <- cbind(di,d0, ">Pi"=tail(subset(tmp,tmp$pvalue<0.05 & tmp$alpha>0),1)[1,1],sign_pos=sum(tmp$pvalue<0.05 & tmp$alpha>0),">Pi"=head(subset(tmp,tmp$pvalue<0.05 & tmp$alpha<0),1)[1,1], sign_neg=sum(tmp$pvalue<0.05 & tmp$alpha<0))

        output[['subset']] <- tmpSubset   
    }

    
    p <- ggplot(tmp, aes(x = Pi)) + geom_point(aes(y = alpha,  colour = pvalue < 0.05)) +scale_y_continuous(breaks = pretty(output$alpha, n = 5)) +
    scaleColourPublication(name="P-value") + labs(y = expression(italic(α)),x = "Pi") + theme(legend.position = c(0.8, 0.9))  + themePublication() + ggtitle(paste("Di: ", di, " - D0: ", d0))
    
    output[['results']] <- tmp
    output[['plot']] <- p

    return(output)
}