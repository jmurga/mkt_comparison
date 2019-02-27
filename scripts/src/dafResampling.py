import pandas as pd
def dafWithResampling(id,data,resamplingValue,type):

	columns = ['id','rawDerivedAllele','div','type']
	output = pd.DataFrame(columns=columns)

	if(data.shape[0] < resamplingValue):
		output = pd.DataFrame({'id':id,'rawDerivedAllele':0,'div':0,'type':type},index=['0'])
		div = output.groupby(['id','type'])['div'].sum().reset_index()
		div = div[['id','div','type']]

		daf = output[['id','rawDerivedAllele','type']]
		return(daf,div)
	else:
		# Delete reference sequence
		ref = data.iloc[[0]]
		outgroup = data.iloc[[-1]]
		data = data.iloc[1:len(data)-1]
		# print(len(data.columns.tolist()))
		
		for j in data.columns.tolist():
			# print(j)
			div = 0
			af = 0
			
			if(data[[j]][data[[j]] != 'N'].dropna().shape[0] < 160):
				continue
			else:
				#Sampling
				tmp = data[[j]][data[[j]]!='N'].dropna().sample(160,replace=False)
				# Merging outgroup
				tmp = pd.concat([tmp,outgroup[[j]]])
				pos = pd.concat([ref[[j]],tmp]).reset_index(drop=True)
				
				if(pos.loc[len(pos)-1,j]=='N' or pos.loc[len(pos)-1,j]=='-'):
					continue
				elif((pos.loc[len(pos)-1,j] != pos.loc[0,j]) & (len(pos.loc[1:len(pos)-2,j].unique())==1)): 
					div = 1
					
				else:

					AA = pos.loc[len(pos)-1,j]
					AN = 160
					AC = pos.loc[1:len(pos)-2,j].value_counts()

					if(AA not in AC.index):
						af=0
					else:
						AC = AC[AC!=AC[AA]]
						if(len(AC) == 0):
							af=0
						else:
							af=AC[0]/AN
			tmp = pd.DataFrame({'id':id,'rawDerivedAllele':af,'div':div,'type':type},index=[id])
			tmp = tmp.reset_index(drop=True)
			output = pd.concat([output,tmp])


		div = output.groupby(['id','type'])['div'].sum().reset_index()
		div = div[['id','div','type']]
		if(type == '4fold'):
			div.columns = ['id','d0','type']
		else:
			div.columns = ['id','di','type']

		daf = output[['id','rawDerivedAllele','type']][output['rawDerivedAllele']!=0]

	return(daf,div)
