daf <- data.frame("Non-synonymous"=c(7,2,1,0,0,0,0,1,0,15),"Synonymous"=c(6,5,2,1,0,1,0,1,1,8))
daf<-melt(daf)
daf$DAF <- rep(seq(0.1,1,0.1),2)
 
p<-ggplot(daf, aes(x = DAF, y = value, fill = variable)) + geom_col(position = "dodge") + theme_bw() + scale_fill_manual(values=c(paletteSanMiguel[3],"gray"),breaks = c("Non.synonymous", "Synonymous"),  labels = c("Non-synonymous", "Synonymous") ) + labs(x="Derived allele frequency", y="Number of sites", fill="Type") + scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + geom_text(aes(label = value), vjust = -0.5, position = position_dodge(width = 0.1)) + ylim(0,20) +theme(legend.position="bottom")

 

ggsave(p,filename='daf.svg',width=10)
