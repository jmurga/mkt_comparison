

d = df[(df_variable !=  'emkt2') & (df_variable !=  'fww2') &  (df_variable !=  'imp2') & (df_variable !=  'asymp2')]


d_variable = d_variable_cat_rename_categories({'std':'MKT','emkt1': 'eMKT 0_05','emkt3': 'eMKT 0_25','fww1': 'fwwMKT 0_05','fww3': 'fwwMKT 0_25','imp1': 'impMKT 0_05','imp3': 'impMKT 0_25','asymp1':'aMKT','grapes':'Grapes','trueAlpha': 'True α'})

p = ggplot(d,aes(x='variable',y='value',fill='variable')) + geom_boxplot(show_legend=False) + facet_wrap('~simulation') + theme_bw() + theme(strip_text=element_text(size=14),axis_text_x = element_text(size=12,angle=90),axis_text_y = element_text(size=12),axis_title=element_text(size=14)) + scale_fill_cmap_d() + ylab("α") + xlab("")  

ggsave(p,filename='/home/jmurga/weakly_pdf',dpi=300) 

#######################

s1 = pd.read_csv("/home/jmurga/mkt/201902/rawData/simulations/base_weakly/sfs.tsv",sep='\t')
s1['pi_nopos'] = s1.pi - s1.pw
s1 = cumulativeSfs(s1.to_numpy())
d1= pd.read_csv("/home/jmurga/mkt/201902/rawData/simulations/base_weakly/div.tsv",sep='\t').to_numpy().flatten() 

a1 = amkt(s1,d1)
b1 = amkt(s1[:,[0,4,2]],d1)

df1 = pd.DataFrame({'f':np.arange(1,40)/40,'All alleles':a1[1],'Neutral + deleterious':b1[1],'simulation':'weakly'}) 


s2 = pd.read_csv("/home/jmurga/mkt/201902/rawData/simulations/base_weakly_bgs/sfs.tsv",sep='\t')
s2['pi_nopos'] = s2.pi - s2.pw
s2 = cumulativeSfs(s2.to_numpy())
d2 = pd.read_csv("/home/jmurga/mkt/201902/rawData/simulations/base_weakly_bgs/div.tsv",sep='\t').to_numpy().flatten() 


a2 = amkt(s2,d2)
b2 = amkt(s2[:,[0,4,2]],d2)

df2 = pd.DataFrame({'f':np.arange(1,40)/40,'All alleles':a2[1],'Neutral + deleterious':b2[1],'simulation':'weakly_bgs'}) 

d = pd.melt(pd.concat([df1,df2]),id_vars=['f','simulation']) 
 

p = ggplot(d,aes(x='f',y='value',color='variable')) + geom_line() + facet_wrap('~simulation') + scale_color_manual(values=['black','#ab2710'],name = "SFS", labels = ["All alleles","Neutral + deleterious"]) + theme_bw() + ylab("α") +  xlab('Derived Allele Count') + theme(legend_position="bottom",legend_text=element_text(size=14),plot_title=element_text(hjust=0.5,face='bold'))

ggsave(p,filename='/home/jmurga/weakly_alpha.pdf',dpi=300) 