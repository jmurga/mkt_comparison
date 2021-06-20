using Analytical, CSV, DataFrames, OrderedCollections, PyCall


function mktByGene(data,population,cutoff=0.025)

	p = pyimport("pyAmkt")

	df = CSV.read(data,DataFrame)
	geneList = df.id |> unique;

	output = []
	@showprogress for iter in geneList
		try 
			sfs,d,m,a = parseBinnedSfs(data=data,population=population,geneList=[iter],cumulative=false);
		catch
			continue
		end

		std =  standardMK(sfs=sfs,divergence=d);
		
		checkDaf = sfs[sfs[:,1] .> cutoff,:]

		if(sum(checkDaf[:,2]) == 0 || sum(checkDaf[:,3]) == 0)
			fww = OrderedDict{String,Float64}("alpha"=>NaN,"pvalue"=>NaN)
			imp = OrderedDict{String,Float64}("alpha"=>NaN,"pvalue"=>NaN)
		else
			fww = fwwMK(sfs=sfs,divergence=d,m=m,cutoff=cutoff);
			imp = impMK(sfs=sfs,divergence=d,m=m,cutoff=cutoff);
		end

		cSfs = Analytical.cumulativeSfs(sfs)
		asymp = p.amkt(sfs,d)[1]
		
		alphas = [std["alpha"],fww["alpha"],imp["alpha"],asymp["alpha"]]
		pvalues = [std["pvalue"],fww["pvalue"],imp["pvalue"], 0]  
		tmp = DataFrame([fill(iter,4),fill(population,4),alphas,pvalues,["standardMK","fwwMK","impMK","aMK"]],[:id,:pop,:alpha,:pvalue,:test])
		push!(output,tmp)
	end

	df = reduce(vcat,output);

	alphaTable = []
	for d in groupby(df, :test)
		tmp = filter(row -> ! isnan(row.alpha), d)[!,[:alpha,:pvalue,:id]]
		tmpPositive = tmp[(tmp.alpha .> 0) .& (tmp.alpha .<1 ) .& (tmp.pvalue .< 0.05),:]
		tmpNegative = tmp[(tmp.alpha .< 0) .& (tmp.pvalue .< 0.05),:]

		msd = mean_and_var(tmp[!,:alpha])
		analyzable = ["analyzable" size(tmp,1) msd[1] msd[2] unique(d.test) [tmp.id]]
		msdP = mean_and_var(tmpPositive[!,:alpha])
		positive = ["positive" size(tmpPositive,1) msdP[1] msdP[2] unique(d.test) [tmpPositive.id]]
		msdN = mean_and_var(tmpNegative[!,:alpha])
		negative = ["negative" size(tmpNegative,1) msdN[1] msdN[2] unique(d.test) [tmpNegative.id]]

		push!(alphaTable,vcat(analyzable,positive,negative))
	end

	result = DataFrame(vcat(alphaTable...),[:analysis,:n,:mean,:sd,:test,:id])

	return(output)
end

function samplingGenes(geneList,replicas,sample)

	output = StatsBase.sample.(fill(geneList,replicas),fill(sample,replicas),replace=true)

	return(output)
end

function binByRecombRate(df,geneList,bins)

	grp = string.(cut(df.cM_Mb,5,labels=[string(i) for i in 1:bins]))
	labels = string.(cut(df.cM_Mb,bins))

	insertcols!(df,1,:grp => grp,:lbls=>labels)
	d = groupby(df,:grp)
	output = [unique(Array(i.id)) for i in d[1:end-1]]

	return(output,labels)
end

function sampleAnalysis(data,geneList,population,sampling,replicas,sample,recomb=true)

	if recomb 
		sampling = samplingGenes(geneList,replicas,size(geneList,1))
	else
		sampling = samplingGenes(geneList,replicas,sample)	
	end

	output = []
	@showprogress for i in sampling
		try 
			sfs,d,m,a = parseBinnedSfs(data=data,population=population,geneList=i,cumulative=false);
		catch
			continue
		end

		if (sum(sfs[:,3]) == 0 || d[2]  == 0 || sum(sfs[:,2]) == 0 || d[1] == 0 || m[1]  == 0 || m[2] == 0 || m[2] < sum(sfs[:,2:3]))

			tmp = DataFrame([NaN size(i,1)],[:aMK,:bins])
		else
			
			asymp = p.amkt(sfs,d)[1]
			tmp = DataFrame([asymp["alpha"] size(i,1)],[:aMK,:bins])
		end
		push!(output,tmp)
	end

	df = vcat(output...)

	output = do.call(rbind.data.frame,tmp)
	output[['bin']] = bins
	colnames(output) <- c('alphaStandard','alphaEmkt1','alphaEmkt2','alphaFww1','alphaFww2','alphaAsymptotic1','bins')
	return(output)
end

function binnedDataToDofe()

end

vennPlot <- function(a,b,c,population,palette){
	require('VennDiagram')

	venn.diagram(
		x              = list(a,b,c),
		category.names = c(paste0('Standard (',length(a),')') , paste0('FWW (',length(b),')'), paste0('imputedMKT (',length(c),')')),
		filename       = paste0('/home/jmurga/mktComparison/results/alphaTables/',population,'Diagram.png'),
		output         = TRUE,
		# Output features
		imagetype      = 'png' ,
		height         = 800 , 
		width          = 800 , 
		resolution     = 300,
		compression    = 'lzw',

		# Circles
		lwd  = 2,
		lty  = 'blank',
		fill = palette,

		# Numbers
		cex = .6,
		fontface = 'bold',
		fontfamily = 'sans',
		
		# Set names
		main = paste0(population,' positive selected genes'),
		cat.cex = 0.6,
		cat.fontface = 'bold',
		cat.default.pos = 'outer',
		cat.pos = c(-27, 27, 135),
		cat.dist = c(0.055, 0.055, 0.1),
		cat.col = palette,
		cat.fontfamily = 'sans',
		rotation = 1
	)
}
