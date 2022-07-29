# from numba import njit
import numpy as np
from scipy import optimize
from scipy import stats
import pandas as pd 
from fisher import pvalue
import sys

# DAF SHOULD BE A NUMPY ARRAY CONTAINING FREQ, PI AND P0 IN THAT ORDER
# DIV SHOULD BE A NUMPY ARRAY CONTAINING DI, D0, MI, M0
# CUMULATIVESFS. DESCRIBED AT HTTPS://STATIC-CONTENT.SPRINGER.COM/ESM/ART%3A10.1038%2FS41559-019-0890-6/MEDIAOBJECTS/41559_2019_890_MOESM1_ESM.PDF USING ALL THE INFORMATION ABOVE THE FREQUENCY. THE FIRST CATEGORY AT SFS INCLUDE PI

def amkt(daf, div, xlow=0, xhigh=1):
	output = {}

	dRatio = float(div[1] / div[0])
	# Compute alpha values and trim
	alpha = 1 - dRatio * (daf[:,1] / daf[:,2])
	trim = ((daf[:,0] >= xlow) & (daf[:,0] <= xhigh))

	# Two-step model fit:
	# First bounded fit:
	try:
		popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim],bounds=([-1, -1, 1], [1, 1, 10]))
		# print('fit initial')
	except:
		# print('could not fit initial')
		popt = None
		pcov = None


	# Second fit using initially guessed values or unbounded fit:
	try:
		popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim], p0=popt,method='lm')
		# print('Fit: lm')
	except:
		try:
			popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim], p0=popt,  method='trf')
			# print('Fit: trf')
		except:
			try:
				popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim], p0=popt, method='dogbox')
				# print('Fit: dogbox')
			except:
				popt=None

	if popt is None:
		output = {'a': np.nan,'b':np.nan,'c': np.nan,'alpha':np.nan,'ciLow':np.nan,'ciHigh':np.nan}

	else:
		output['a'] = popt[0]
		output['b'] = popt[1]
		output['c'] = popt[2]

		# alpha for predicted model
		output['alpha'] = exp_model(1.0, output['a'], output['b'], output['c'])
		
		# Compute confidence intervals based on simulated data (MC-SOERP)
		vcov = pd.concat([pd.DataFrame([0] * 4).transpose(),
						  pd.concat([pd.DataFrame([0] * 4), pd.DataFrame(pcov)], axis=1, ignore_index=True)],
						 axis=0, ignore_index=True)
		vcov = vcov.iloc[0:4, :].values

		try:
			simpars = np.random.multivariate_normal(mean=[1.0, output['a'], output['b'], output['c']], cov=vcov, size=10000, check_valid='raise')  # same as R implementation
			output['ciLow'], output['ciHigh'] = np.quantile([exp_model(x[0], x[1], x[2], x[3]) for x in simpars], [0.025, 0.975])
		except:
			output['ciLow'], output['ciHigh']= np.nan, np.nan

	return output,alpha[trim]

def exp_model(f_trimmed, a, b, c):
	return a + b * np.exp(-c * f_trimmed)

def cumulativeSfs(x):

	f       = x[:,0]
	sfsTemp = x[:,1:]
	
	out       = np.empty_like(x)
	out[0,0]  = f[0]
	out[0,1:] = np.sum(sfsTemp,axis=0)

	for i in range(1,out.shape[0]):

		app = out[i-1,1:] - sfsTemp[i-1,:]

		if np.sum(app) > 0.0:
			out[i,0]  = f[i]
			out[i,1:] = app
		else:
			out[i,0]  = f[i]
			out[i,1:] = np.zeros(app.shape[0])

	return out

def reduceSfs(x,bins):

	bins = (bins*2) - 1 
	f = x[:,0]
	sfs = x[:,1:]
	
	b = np.arange(0,1,1/bins)
	inds = np.digitize(f,b,right=True)
	out  = np.zeros((bins,x.shape[1]))
	out[:,0]  = np.unique(inds)
	
	sfsGrouped = np.hstack([np.reshape(inds,(inds.shape[0],1)),sfs])
	for i in np.unique(inds):
		out[out[:,0]==i,1:] = np.sum(sfsGrouped[sfsGrouped[:,0] == i,1:],axis=0) 
		
	out[:,0] = b
	return(out)

def impMK(sfs, divergence,l,h=None,m=None):

	output = {}

	pn = np.sum(sfs[:,1])
	ps = np.sum(sfs[:,2])
	dn = divergence[0]
	ds = divergence[1]


	toFix = 0
	deleterious = 0
	pnHigh = 0
	psHigh = 0
	### Estimating slightly deleterious with pn/ps ratio
	fltLow = (sfs[:, 0] <= l)
	pnLow   = sfs[fltLow][:,1].sum()
	psLow   = sfs[fltLow][:,2].sum()

	if (h is None):
		fltInter = (sfs[:, 0] >= l) & (sfs[:, 0] <= 1)
		pnInter = sfs[fltInter][:, 1].sum()
		psInter = sfs[fltInter][:, 2].sum()
	else:
		fltInter = (sfs[:, 0] >= l) & (sfs[:, 0] < h)
		pnInter = sfs[fltInter][:, 1].sum()
		psInter = sfs[fltInter][:, 2].sum()

		fltHigh = (sfs[:, 0] >= h)
		pnHigh  = sfs[fltHigh][:, 1].sum()
		psHigh  = sfs[fltHigh][:, 2].sum()

		toFix = pnHigh - (pnInter*psHigh/psInter)
		if toFix < 0:
			toFix = 0
		dn = round(dn + abs(toFix),3)


	ratiops       = psLow / psInter
	deleterious   = pnLow - (pnInter * ratiops)
	pnNeutral     = round(pn - deleterious - toFix,3)

	if(deleterious < 0):
		deleterious = 0

	output['alpha'] = round(1 - (((pn - deleterious) / ps) * (ds / dn)),3)
	# output[0] = round(1 - ((pnNeutral/ps) * (ds / dn)),3)

	output['pvalue'] = pvalue(ps, ds, pnNeutral, dn).two_tail
	output['s'] = pn+ps	
	output['deleterious'] = deleterious

	if(m is not None):
		m0 = m[1];mi = m[0]
		# ## Estimation of b: weakly deleterious
		output['b'] = (deleterious / ps) * (m0 / mi)

		## Estimation of f: neutral sites
		output['f'] = (m0 * pnNeutral) / (mi * ps)

		## Estimation of d, strongly deleterious sites
		output['d'] = 1 - (output['f'] + output['b'])


		# divergence metrics
		output['Ka']       = dn / mi
		output['Ks']       = ds / m0
		output['omega']    = output['Ka'] / output['Ks']
		output['gamma'] = (pnNeutral/ps - dn/ds) * m0/mi

	return output

def FWW(sfs, divergence,m, cutoff=0.15):

	output = {}

	pn = np.sum(sfs[:,1])
	ps = np.sum(sfs[:,2])
	dn = divergence[0]
	ds = divergence[1]
	mn = m[0]
	ms = m[1]

	### Estimating alpha with pn/ps ratio
	pnGreater = sfs[sfs[:,0] > cutoff,1].sum()
	psGreater = sfs[sfs[:,0] > cutoff,2].sum()

	output['alpha'] = 1 - (pnGreater / psGreater * (ds / dn))
	output['pvalue'] = pvalue(psGreater, ds, pnGreater, dn).two_tail
	output['s'] = pn+ps
	output['deleterious'] = sfs[sfs[:,0] <= cutoff,1:].sum()
	return output

def eMKT(sfs, divergence, m,cutoff=0.15):

	output = {}

	daf_below_cutoff = np.sum(sfs[sfs[:,0] <= cutoff,1:],axis=0) 
	daf_above_cutoff = np.sum(sfs[sfs[:,0] > cutoff,1:],axis=0)
	Pi = sum(sfs[:,1])
	P0 = sum(sfs[:,2]) 
 
	## Estimate fractions
	f_neutral = daf_below_cutoff[1]/P0
	Pi_neutral_below_cutoff =  Pi * f_neutral
	Pi_wd = daf_below_cutoff[0] - Pi_neutral_below_cutoff
	Pi_neutral =  round(Pi_neutral_below_cutoff + daf_above_cutoff[0])
	
	## Estimation of alpha
	output["alpha"] = 1-((Pi_neutral/P0)*(divergence[1]/divergence[0]))
	
	## Estimation of b: weakly deleterious
	b = (Pi_wd/P0)*(m[1]/m[0])
	
	## Estimation of f: neutral sites
	f = (m[1]*Pi_neutral)/(m[0]*P0)
	
	## Estimation of d, strongly deleterious sites
	d = 1-(f+b)
	
	## Fisher exact test p-value from the MKT
	output['pvalue'] = pvalue(P0, divergence[1], Pi_neutral, divergence[0]).two_tail

	output['s'] = Pi + P0
	output['deleterious'] = Pi_wd

	return(output) 

def standardMK(sfs,divergence,m):

	output = {}

	pn = np.sum(sfs[:,1])
	ps = np.sum(sfs[:,2])
	dn = divergence[0]
	ds = divergence[1]
	
	
	output["alpha"] = round(1 - ((pn/ps) * (ds / dn)),3)
	#  method = :mnnlike same results R, python two.sides
	output["pvalue"] = pvalue(ps,pn,ds,dn).two_tail

	mn = m[0]; ms = m[1]

	ka      = dn / mn
	ks       = ds / ms
	output["omega"]    = ka/ ks


	# Omega A and Omega D
	output["omegaA"] = output["omega"] * output["alpha"]
	output["omegaD"] = output["omega"] - output["omegaA"]	    

	return output
