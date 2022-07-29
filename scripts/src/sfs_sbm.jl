
using DataFrames,RCall,MKtest
out = DataFrame[]
amk = DataFrame[]
for i in [5,10,50,100,250,500]

	adap = MKtest.parameters(n=10,gL=i,al_low=0.2,al_tot=0.4)

	x,y = MKtest.analytical_alpha(adap)
	a1  = MKtest.aMK(x)[1]
	a2  = MKtest.aMK(y)[1]

	tmp = DataFrame(:f=>1:19,:All_alleles=>x,:a_x_nopos=>y,:s=>i)
	tmp_amk = DataFrame(:s=>i,:amk_a_x=>a1,:amk_a_x_nopos=>a2)

	push!(out,tmp)
	push!(out,tmp)

end

@rput out

R"""
df = rbindlist(out)
names(df)[2] = "All alleles"
d = melt(df,id.vars=c('f','s'))
d[variable=="a_x"]$variable = "All alleles"
d[variable=="a_x_nopos"]$variable = "Neutral + deleterious alleles"

d$variable = as.factor(d$variable)
d$s = factor(d$s,labels=paste0("2Nes=",unique(d$s)))


p = ggplot(d,aes(x=f,y=value,color=variable)) + geom_line() + facet_wrap(~s) + scale_x_log10() + theme_bw() + scale_color_manual(values=c("black","#ab2710")) +  labs(color="Derived Allele Frequency") + ylab(expression(alpha[x])) + theme(legend.position="bottom")

p

ggsave(p,filename="/home/jmurga/mkt/201902/results/pos.svg")
"""
