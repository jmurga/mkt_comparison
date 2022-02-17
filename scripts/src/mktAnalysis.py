import pandas as pd
import numpy as np
import pickle as pkl
import warnings
from tqdm import tqdm
from scipy import optimize
from fisher import pvalue
from multiprocessing import Pool
from plotnine import *
import subprocess
import glob
import sys
import os
sys.path.insert(0, '/home/jmurga/mkt/201902/scripts/src/')
from mkt import *


# DAF SHOULD BE A NUMPY ARRAY CONTAINING FREQ, PI AND P0 IN THAT ORDER
# DIV SHOULD BE A NUMPY ARRAY CONTAINING DI, D0, MI, M0
# CUMULATIVESFS. DESCRIBED AT HTTPS://STATIC-CONTENT.SPRINGER.COM/ESM/ART%3A10.1038%2FS41559-019-0890-6/MEDIAOBJECTS/41559_2019_890_MOESM1_ESM.PDF USING ALL THE INFORMATION ABOVE THE FREQUENCY. THE FIRST CATEGORY AT SFS INCLUDE PI

def grapesGeneListToHtcondor(path,population,chunk):

	n = np.array(glob.glob(path + "/*.dofe"))

	def getGeneName(ng):
		return(ng.split('/')[-1].split('.')[0])

	vfunc = np.vectorize(getGeneName)

	genes = pd.DataFrame(vfunc(n))
	n = 1999  #chunk row size
	list_df = [genes[i:i+n] for i in range(0,genes.shape[0],n)]
	for i in range(1,len(df) + 1):
		list_df[i].to_csv(PATH + "/rawData/dofe/"+ population +"/list"+str(i)+".txt",index=False)
		
def sfsToDofe(sfs,d,m,output,gene=None):

	if gene:
		name = gene
	else:
		name = 'dataset_bins_' + str(sfs.shape[0])

	data = pd.DataFrame(np.hstack([name,sfs.shape[0] + 1, m[0],sfs[:,1],m[1],sfs[:,2],m[0],d[0],m[1],d[1]])).T

	f = open(output,"w")
	f.write(" \n#unfolded\n")
	f.close()
	data.to_csv(output,sep='\t',header=None,index=False,mode='a')

def parseBinnedSfs(data,population,gList=None):

	tmp  = data[data['pop'] == population]
	
	if(gList is not None):
		tmp = tmp[tmp.id.isin(gList)]

	ds   = sum(tmp.d0)
	dn   = sum(tmp.di)

	ms   = sum(tmp.m0)
	mn   = sum(tmp.mi)

	sfsPs = np.sum(tmp.daf4f.apply(lambda row: np.array(list(map(int,row.split(';'))))))
	sfsPn = np.sum(tmp.daf0f.apply(lambda row: np.array(list(map(int,row.split(';'))))))

	f = np.arange(1,sfsPs.shape[0]+1)/ sfsPs.shape[0]
	sfs = np.stack([f,sfsPn,sfsPs],axis=1)

	if len(gList) > 1:
		g = gList
	else:
		g = gList[0]

	return(sfs,np.array([dn,ds]),np.array([mn,ms]),g)

def runGrapes(inp,out):
	subprocess.run(["grapes", "-in",inp,"-out",out,"-model","GammaZero"])

def dataEstimates(sfs,d,m,geneId,population,dofe):

	# Manual cleaning to input names easly in HTcondor
	# Gene-by-gene mantaining the ID
	# Sampling changing to Iter_concat.dofe. If sampling == 1 then str and change a posteriori
	if type(geneId) == str:
		if dofe.split('/')[-1] != '':
			dofeFile = dofe + '_concat.dofe'
			dataset = dofe.split('/')[-1] + '_concat'

		else:
			dofeFile = dofe + geneId + '.dofe'
			dataset = geneId
	else:
		dofeFile = dofe + '_concat.dofe'
		dataset = dofe.split('/')[-1] + '_concat'


	sfsToDofe(sfs = sfs,d = d,m = m,output = dofeFile,gene=dataset)

	# Execute each mkt test through Analytical.jl
	## StandardMKT
	std =  standardMK(sfs=sfs,divergence=d,m=m);

	checkDaf = sfs[sfs[:,0] > 0.15,:]
	if((checkDaf[:,1].sum() == 0) | (checkDaf[:,2].sum() == 0)):
		fww1 = {"alpha":np.nan,"pvalue":np.nan,'s': np.nan, 'deleterious': np.nan}
		imp1 = {"alpha":np.nan,"pvalue":np.nan,'s': np.nan, 'deleterious': np.nan}

	else:
		fww1 = FWW(sfs=sfs,divergence=d,m=m,cutoff=0.15);
		imp1 = impMK(sfs=sfs,divergence=d,m=m,l=0.15);


	if(checkDaf[:,2].sum() == 0):
		emkt1 = {"alpha":np.nan,"pvalue":np.nan,'s': np.nan, 'deleterious': np.nan}
	else:
		emkt1 = eMKT(sfs=sfs,divergence=d,m=m,cutoff=0.15)
	
	checkDaf = sfs[sfs[:,0] > 0.25,:]
	if((np.sum(checkDaf[:,1]) == 0) | (np.sum(checkDaf[:,2]) == 0)):
		fww2 = {"alpha":np.nan,"pvalue":np.nan,'s': np.nan, 'deleterious': np.nan}
		imp2 = {"alpha":np.nan,"pvalue":np.nan,'s': np.nan, 'deleterious': np.nan}
	else:
		fww2 = FWW(sfs=sfs,divergence=d,m=m,cutoff=0.25);
		imp2 = impMK(sfs=sfs,divergence=d,m=m,l=0.25);

	if(checkDaf[:,2].sum() == 0):
		emkt2 = {"alpha":np.nan,"pvalue":np.nan,'s': np.nan, 'deleterious': np.nan}
	else:
		emkt2 = eMKT(sfs=sfs,divergence=d,m=m,cutoff=0.25)


	if (sfs == 0).any():
		asymp1 = amkt(sfs,d)[0]
		asymp2 = amkt(sfs,d,0.1,0.9)[0]
	else:
		cSfs = cumulativeSfs(sfs)
		asymp1 = amkt(cSfs,d)[0]
		asymp2 = amkt(cSfs,d,0.1,0.9)[0]

	pvalues = np.array([std["pvalue"],emkt1["pvalue"],emkt2["pvalue"],fww1["pvalue"],fww2["pvalue"],imp1["pvalue"],imp2["pvalue"],0,0])
	alphas = np.array([std["alpha"],emkt1["alpha"],emkt2["alpha"],fww1["alpha"],fww2["alpha"],imp1["alpha"],imp2["alpha"],asymp1["alpha"],asymp2["alpha"]])
	n = alphas.shape[0]

	tmp = pd.DataFrame([np.repeat(dataset,n),np.repeat(population,n),alphas,pvalues,np.array(["std","emkt1","emkt2","fww1","fww2","imp1","imp2","asymp1","asymp2"])]).T
	tmp.columns= ["id","pop","alpha","pvalue","test"]
	deleterious = np.array([0,emkt1["deleterious"],emkt2["deleterious"],fww1["deleterious"],fww2["deleterious"],imp1["deleterious"],imp2["deleterious"]])
	s = np.array([sfs[:,1].sum(),emkt1["s"],emkt2["s"],fww1["s"],fww2["s"],imp1["s"],imp2["s"]])
	tmp_deleterious = pd.DataFrame([np.repeat(dataset,6),np.repeat(population,6),s,deleterious,np.array(["std","emkt1","emkt2","fww1","fww2","imp1","imp2"])]).T
	tmp_deleterious.columns= ["id","pop","s","deleterious","test"]


	return(tmp,tmp_deleterious)

def mktByGene(df,population,nthreads,dofe,grapes=False):

	outDofe = dofe + "/" + population.lower() + "/"
	os.makedirs(outDofe,exist_ok=True)

	tmp = df[df['pop'] == population]
	n = tmp.shape[0]
	geneList = tmp.id.unique()
	geneList = [[i] for i in geneList.tolist()]

	pool = Pool(processes = nthreads)
	lSfs,lDiv,lM,genes = zip(*pool.starmap(parseBinnedSfs,zip([tmp]*n,[population]*n,geneList)))
	pool.terminate()

	pool = Pool(processes = nthreads)
	tmpAlphas,tmp_deleterious = zip(*pool.starmap(dataEstimates,zip(lSfs,lDiv,lM,genes,[population]*n,[outDofe]*n)))
	pool.terminate()	

	if grapes:
		cmd = "ls " + outDofe + " | cut -d'.' -f1 | parallel -j" + str(nthreads) + " 'timeout 10m grapes -in {}.dofe -out {}.out -model GammaZero'"
		subprocess.run(cmd,shell=True)

	dfAlpha = pd.concat(tmpAlphas);
	dfAlpha[dfAlpha.alpha > 1] = np.nan
	dfAlpha.alpha = dfAlpha.alpha.astype(float)
	dfAlpha.pvalue = dfAlpha.pvalue.astype(float)
	flt = dfAlpha

	analyzable = flt.groupby(['test','pop'], as_index=False).agg({'alpha':['count','mean','std'],'id':','.join}).reset_index(drop=True)
	analyzable.columns = analyzable.columns.droplevel(0) 
	analyzable.columns = ['test','pop','count', 'mean', 'std','ids']

	# positive = flt[(flt.alpha > 0)].groupby( ['test','pop'], as_index=False).agg({'alpha':['count','mean','std'],'id':','.join}).reset_index(drop=True)
	positive = flt[(flt.alpha > 0) & (flt.pvalue <= 0.05)].groupby( ['test','pop'], as_index=False).agg({'alpha':['count','mean','std'],'id':','.join}).reset_index(drop=True)
	positive.columns = positive.columns.droplevel(0) 
	positive.columns = ['test','pop','count', 'mean', 'std', 'ids']

	negative = flt[(flt.alpha < 0) & (flt.pvalue <= 0.05)].groupby( ['test','pop'], as_index=False).agg({'alpha':['count','mean','std'],'id':','.join}).reset_index(drop=True)
	negative.columns = negative.columns.droplevel(0) 
	negative.columns = ['test','pop','count', 'mean', 'std','ids']

	analyzable['type'] = 'analyzable'
	positive['type'] = 'positive'
	negative['type'] = 'negative'

	positiveRaw = flt[(flt.alpha > 0)]

	df_deleterious = pd.concat(tmp_deleterious).infer_objects()

	return(flt,pd.concat([analyzable,positive,negative]),df_deleterious)

def samplingGenes(geneList,replicas,sample):

	output = np.array([np.random.choice(geneList,size=sample,replace=True) for i in range(0,replicas)])
	return(output)

def binByRecombRate(dfRecomb,replicas,selectBin,bins):

	dfRecomb['grp'] = pd.qcut(dfRecomb.cM_Mb,bins,duplicates='drop')
	# tmpRecomb = [dfRecomb[dfRecomb.grp == i].id.unique().tolist() for i in dfRecomb.grp.cat.categories]

	tmpRecomb = dfRecomb[dfRecomb.grp == dfRecomb.grp.cat.categories[selectBin]].id.unique()

	gList = [np.random.choice(tmpRecomb,900,replace=True) for i in range(0,replicas)]
	
	grp = str(dfRecomb.grp.cat.categories[selectBin])

	return(gList,grp)

def sampleAnalysis(df,population,sample,replicas,dofe,nthreads,grapes=False,recombBin=None,bins=None):

	if bins is not None:
		outDofe = dofe + "/" + population.lower() + "/" + str(recombBin)
		os.makedirs(outDofe, exist_ok=True)
		outDofeList = [ outDofe + "/" + str(i) for i in range(1,replicas+1)]

		tmp = df[df['pop'] == population]
		sampling, grp = binByRecombRate(tmp,replicas,recombBin,bins)
	else:

		tmp = df[df['pop'] == population]
		allGenes = tmp.id.unique()
		np.random.seed(1331272)
		geneList = np.random.choice(allGenes,3500)

		outDofe = dofe + "/" + population.lower() + "/" + str(sample)
		os.makedirs(outDofe, exist_ok=True)
		outDofeList = [ outDofe + "/" + str(i) for i in range(1,replicas+1)]
		sampling = samplingGenes(geneList,replicas,sample)

	pool = Pool(processes = nthreads)
	lSfs,lDiv,lM,genes = zip(*pool.starmap(parseBinnedSfs,zip([tmp]*replicas,[population]*replicas,sampling)))
	pool.terminate()

	pool = Pool(processes = nthreads)
	tmpAlphas = pool.starmap(dataEstimates,zip(lSfs,lDiv,lM,genes,[population]*replicas,outDofeList))
	pool.terminate()

	if bins is not None:

		dfAlpha = pd.concat(tmpAlphas)
		dfAlpha[dfAlpha.alpha > 1] = np.nan
		dfAlpha.replace([np.inf, -np.inf], np.nan,inplace=True)
		dfAlpha.alpha = dfAlpha.alpha.astype(float)
		dfAlpha.pvalue = dfAlpha.pvalue.astype(float)
		flt = dfAlpha.dropna()

		analyzable = flt.groupby(['test','pop'], as_index=False).agg({'alpha':['count','mean','std']}).reset_index(drop=True)
		analyzable.columns = analyzable.columns.droplevel(0) 
		analyzable.columns = ['test','pop','count', 'mean', 'std']
		dfAlpha['bins'] = grp; analyzable['bins'] = grp
	else:	
		dfAlpha = pd.concat(tmpAlphas)
		dfAlpha[dfAlpha.alpha > 1] = np.nan
		dfAlpha.replace([np.inf, -np.inf], np.nan,inplace=True)
		dfAlpha.alpha = dfAlpha.alpha.astype(float)
		dfAlpha.pvalue = dfAlpha.pvalue.astype(float)
		
		analyzable = dfAlpha.groupby(['test','pop'], as_index=False).agg({'alpha':['count','mean','std']}).reset_index(drop=True)
		analyzable.columns = analyzable.columns.droplevel(0) 
		analyzable.columns = ['test','pop','count', 'mean', 'std']
		dfAlpha['bins'] = str(sample) 
		analyzable['bins'] = str(sample)

	return(dfAlpha,analyzable)

def mktOnsimulatedData(path,model,nthreads):

	sfsFiles = np.sort(glob.glob(path + "/sfs*"))[1:]
	divFiles = np.sort(glob.glob(path + "/div*"))[2:]
	dofeFiles = np.sort(glob.glob(path + "/*" + model))

	pool = Pool(processes = nthreads)
	alphas = pool.starmap(simulationsEstimates,zip(sfsFiles,divFiles,dofeFiles))
	pool.terminate()

	df = pd.DataFrame(np.vstack(alphas),columns=['std','emkt1','emkt2','emkt3','emkt4','fww1','fww2','fww3','fww4','imp1','imp2','imp3','imp4','asymp1','asymp2','grapes','trueAlpha'])
	df = np.round(df,3)
	df['simulation'] = path.split('/')[-1]

	dfPlot = pd.melt(df,id_vars=['simulation'])
	dfPlot.variable = pd.Categorical(dfPlot.variable, categories=dfPlot.variable.unique(), ordered=True)

	p = ggplot(dfPlot,aes(x='simulation',y='value',fill='variable')) + geom_boxplot()  
	std_error = np.mean(abs(df.trueAlpha - df['std']))
	emkt1_error = np.mean(abs(df.trueAlpha - df.emkt1))
	emkt2_error = np.mean(abs(df.trueAlpha - df.emkt2))
	emkt3_error = np.mean(abs(df.trueAlpha - df.emkt3))
	fww1_error = np.mean(abs(df.trueAlpha - df.fww1))
	fww2_error = np.mean(abs(df.trueAlpha - df.fww2))
	fww3_error = np.mean(abs(df.trueAlpha - df.fww3))
	imp1_error = np.mean(abs(df.trueAlpha - df.imp1))
	imp2_error = np.mean(abs(df.trueAlpha - df.imp2))
	imp3_error = np.mean(abs(df.trueAlpha - df.imp3))
	asymp1_error = np.mean(abs(df.trueAlpha - df.asymp1))
	asymp2_error = np.mean(abs(df.trueAlpha - df.asymp2))
	grapes_error = np.mean(abs(df.trueAlpha - df.grapes))


	error = pd.DataFrame([std_error,emkt1_error,emkt2_error,emkt3_error,fww1_error,fww2_error,fww3_error,imp1_error,imp2_error,imp3_error,asymp1_error,asymp2_error,grapes_error],columns=['value'])
	error['test'] = ['std','emkt1','emkt2','emkt3','fww1','fww2','fww3','imp1','imp2','imp3','asymp1','asymp2','grapes']
	error['stat'] = 'error'

	meanSd = df.groupby(['simulation'], as_index=False).agg(['mean','std']).T.reset_index()
	meanSd.columns =  ["test","stat","value"]

	errorMeanSd = pd.concat([error,meanSd])
	return(df,dfPlot,errorMeanSd,p)

def simulationsEstimates(sfile,dfile,outGrapes):

	sfs = pd.read_csv(sfile,sep='\t').to_numpy()
	d   = pd.read_csv(dfile,sep='\t').to_numpy().flatten()
	m   = np.array([21*10**6*0.75,21*10**6*0.25])

	if d.shape[0] > 3:
		trueAlpha = np.round(np.sum(d[2:])/d[0],3)

	else:
		trueAlpha = np.round(d[2]/d[0],3)
	# Execute each mkt test through Analytical.jl
	## StandardMKT
	std =  standardMK(sfs=sfs,divergence=d,m=m);

	cSfs = cumulativeSfs(sfs)

	emkt1 = eMKT(sfs=sfs,divergence=d,m=m,cutoff=0.05)
	emkt2 = eMKT(sfs=sfs,divergence=d,m=m,cutoff=0.15)
	emkt3 = eMKT(sfs=sfs,divergence=d,m=m,cutoff=0.35)
	emkt4 = eMKT(sfs=sfs,divergence=d,m=m,cutoff=0.35)

	fww1  = FWW(sfs=sfs,divergence=d,m=m,cutoff=0.05);
	fww2  = FWW(sfs=sfs,divergence=d,m=m,cutoff=0.15);
	fww3  = FWW(sfs=sfs,divergence=d,m=m,cutoff=0.25);
	fww4  = FWW(sfs=sfs,divergence=d,m=m,cutoff=0.35);

	imp1  = impMK(sfs=sfs,divergence=d,m=m,l=0.05);
	imp2  = impMK(sfs=sfs,divergence=d,m=m,l=0.15);	
	imp3  = impMK(sfs=sfs,divergence=d,m=m,l=0.25);
	imp4  = impMK(sfs=sfs,divergence=d,m=m,l=0.35);


	asymp1 = amkt(cSfs,d)[0]
	asymp2 = amkt(cSfs,d,0.1,0.9)[0]

	# subprocess.run(["grapes", "-in",dofe,"-out",outGrapes,"-model","GammaZero"])
	
	grapes = pd.read_csv(outGrapes)

	alphas = np.array([std["alpha"],emkt1["alpha"],emkt2["alpha"],emkt3["alpha"],emkt4["alpha"],fww1["alpha"],fww2["alpha"],fww3["alpha"],fww4["alpha"],imp1["alpha"],imp2["alpha"],imp3["alpha"],imp4["alpha"],asymp1["alpha"],asymp2["alpha"],grapes.alpha.iloc[1],trueAlpha])
	return alphas

def openFiles(sfsFile,divFile):

	sfs = pd.read_csv(sfsFile,sep='\t').to_numpy()
	div = pd.read_csv(divFile,sep='\t').to_numpy().flatten()

	return(sfs,div)

def parseBootstrapPolDiv(path,N,genes,sample=1,replicas=100,nthreads=8):
	"""
	slr, tupple array with daf and div data by element in list
	"""
	dafFiles = np.sort(glob.glob(path + "/daf/*.tsv.gz"))
	divFiles = np.sort(glob.glob(path + "/div/div*.tsv.gz"))

	idxAll = np.arange(0,len(dafFiles)-1)
	bn = (N*2) - 1
	rsample = int(sample * len(dafFiles))

	pool = Pool(processes = nthreads)
	lSfs,ldiv = zip(*pool.starmap(openFiles,list(zip(dafFiles,divFiles))))
	pool.terminate()

	bn = (N*2)-1
	sfs = np.sum(lSfs,axis=0)[:,1:]
	f = np.reshape(np.arange(1,bn+1)/bn,(bn,1))
	sfs = np.hstack((f,sfs))
	# sCumu = cumulativeSfs(sfs)
	d = np.sum(ldiv,axis=0)

	if sfs.shape[1] > 3:
		sC = ['f','pi','p0','pw']
		dC = ['di','d0','dw','ds']
	else:
		sC = ['f','pi','p0']
		dC = ['di','d0','ds']

	fName      = path + path.split("/")[-1]
	sfs        = pd.DataFrame(np.round(sfs,5),columns=sC)
	d          = pd.DataFrame(d).T
	d.columns  = dC

	sfs.to_csv(path + "/sfs.tsv", header=True, index=False, sep="\t")
	d.to_csv(path + "/div.tsv", header=True, index=False, sep="\t")

	for r in tqdm(range(1, replicas+1)):
		idx    = np.sort(np.random.choice(idxAll,rsample,replace=True))
		tmpDaf = [lSfs[i] for i in idx]
		tmpDiv = [ldiv[i] for i in idx]
		m      = np.array([500*3*genes*0.75,500*3*genes*0.25])

		sfs = np.sum(tmpDaf,axis=0)[:,1:]
		f = np.reshape(np.arange(1,bn+1)/bn,(bn,1))
		sfs = np.hstack((f,sfs))
		# sCumu = cumulativeSfs(sfs)
		d = np.sum(tmpDiv,axis=0)
		
		outFile = path + "/" + "simulation" + str(r) + ".dofe"
		sfsToDofe(sfs,d,m,outFile)

		sfs = pd.DataFrame(np.round(sfs,5),columns=sC)
		d = pd.DataFrame(d).T
		d.columns=dC
		sfs.to_csv(path + "/sfs" + str(r) + ".tsv", header=True, index=False, sep="\t")
		d.to_csv(path + "/div" + str(r) + ".tsv", header=True, index=False, sep="\t")

def grapesOutput(path,model,population,nthreads):
	# Assuming extension file == model

	grapes = []
	flist = glob.glob(path + "/*" + model)
	pool = Pool(processes = nthreads)
	n = pool.starmap(readGrapes,zip(flist,[model]*len(flist)))
	pool.terminate()

	df = pd.concat(n)

	tmpGrapes = df[(df.alpha_down < 0.25) & (df.alpha_range !=0)];
	tmp = tmpGrapes.iloc[:,[0,2,6]].values
	tmp = pd.DataFrame({'id':tmp[:,0],'pop':population,'alpha':tmp[:,1].astype(float),'pvalue':np.nan,'test':tmp[:,2]})

	analyzable = tmp.groupby(['test','pop'], as_index=False).agg({'alpha':['count','mean','std'],'id':','.join}).reset_index(drop=True)
	analyzable.columns = analyzable.columns.droplevel(0) 
	analyzable.columns = ['test','pop','count', 'mean', 'std','ids']

	positive = tmp[(tmp.alpha > 0)].groupby( ['test','pop'], as_index=False).agg({'alpha':['count','mean','std'],'id':','.join}).reset_index(drop=True)

	positive.columns = positive.columns.droplevel(0) 
	positive.columns = ['test','pop','count', 'mean', 'std', 'ids']

	negative = tmp[(tmp.alpha < 0)].groupby( ['test','pop'], as_index=False).agg({'alpha':['count','mean','std'],'id':','.join}).reset_index(drop=True)
	negative.columns = negative.columns.droplevel(0) 
	negative.columns = ['test','pop','count', 'mean', 'std','ids']

	analyzable['type'] = 'analyzable'
	positive['type'] = 'positive'
	negative['type'] = 'negative'

	return(df,tmp,pd.concat([analyzable,positive,negative]))

def readGrapes(file,model):
	df = pd.read_csv(file)
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		tmp = df[df.model==model]
		tmp['alpha_up'] = tmp['alpha_up'].astype(float)
		tmp.alpha_down = tmp.alpha_down.astype(float) 
		tmp.omegaA = tmp.omegaA.astype(float) 

	out = tmp.loc[:,['dataset','model','alpha','alpha_down','alpha_up','omegaA']]
	out['alpha_range'] = tmp.alpha_up - tmp.alpha_down
	out['test'] = 'grapes'
	return(out)	

def loadPickleAlphas(path):
	f = glob.glob(path)
	alphas = {}
	binned = {}

	for i in f:
		population = i.split('/')[-1].split('Pooled')[0]
		data = pkl.load(open(i,'rb'))
		a = data[0]
		b = data[1]

		a.bins = a.bins.astype(str)
		b.bins = b.bins.astype(str)        
		
		a.columns = ['geneid', 'pop', 'alpha', 'pvalue', 'test', 'bins']
		a.bins = pd.Categorical(a.bins,categories = a.bins.unique())

		aMelted = pd.melt(a,id_vars=['test','bins','pop'], value_vars='alpha')
		aMelted.test = aMelted.test.astype(str)
		aMelted[['pop']] = aMelted[['pop']].astype(str)
		
		b.bins = pd.Categorical(b.bins,categories = b.bins.unique())
		b.test = b.test.astype(str)
		b['pop'] = b['pop'].astype(str)
		
		alphas[population] = aMelted
		binned[population] = b
	return(alphas,binned)

def grapesBootstrap(bins,replicas,path):

	for i in tqdm(bins):

		o = path + "/" + str(i) + "/bootstrap"
		os.makedirs(o, exist_ok=True)

		df    = pd.read_csv(path + "/" + str(i) + "/429_concat.dofe",sep='\t',skiprows=2,header=None)
		name  = df.iloc[:,0:2].to_numpy()
		m     = df.iloc[:,2:].to_numpy()
		out=[]
		for j in range(1,replicas+1):
			mpois = np.random.poisson(m)
			data = pd.DataFrame(np.hstack([name,mpois]))
			oDofe = o + "/" + str(j) + "_bootstrap.dofe"
			f = open(oDofe,"w")
			f.write(" \n#unfolded\n")
			f.close()
			try:
				out.append(dofeToSfs(data))
			except:
				out.append(np.nan)
			data.to_csv(oDofe,sep='\t',header=None,index=False,mode='a')

def dofeToSfs(df):

	n = df.iloc[:,1].values[0]

	pn = df.iloc[:,3:(n+2)].to_numpy().flatten()
	ps = df.iloc[:,(n+3):(n+n+2)].to_numpy().flatten()
	m = df.iloc[:,n+n+2:].to_numpy().flatten()
	d = m[[1,3]]

	sfs = pd.DataFrame({'f':np.arange(1,pn.shape[0]+1)/pn.shape[0],'pn':pn,'ps':ps}).to_numpy()

	checkDaf = sfs[sfs[:,0] > 0.25,:]
	if((np.sum(checkDaf[:,1]) == 0) | (np.sum(checkDaf[:,2]) == 0)):
		imp1 = {"alpha":np.nan,"pvalue":np.nan}
	else:
		imp1 = impMK(sfs=sfs,divergence=d,m=m,l=0.15);

	return(imp1['alpha'])