import time 
start = time.time()
from pyfaidx import Fasta
file = Fasta('2l.seq')
coordinates = pd.read_csv('test.bed',sep='\t',header=None).iloc[:,[1,2]].values.tolist() 

samples = list(file.keys())

matrix4f = np.empty([len(samples),len(file.get_spliced_seq(samples[0], coordinates))+1],dtype='str')
matrix0f = np.empty([len(samples),len(file.get_spliced_seq(samples[0], coordinates))+1],dtype='str')

# matrix[0] = file.get_spliced_seq(samples[0], coordinates).seq.upper()
refSeq0f,refSeq4f = degenerate(file.get_spliced_seq(samples[0], coordinates).seq.upper()+'C')

deleteIndex = list()
for i in range(1,len(samples),1):
	tmp = file.get_spliced_seq(samples[i],coordinates).seq.upper()
	if(tmp == ('N' * len(tmp))):
		deleteIndex.append(i)
	else:
		matrix4f[i] = list(tmp+'C')
		matrix0f[i] = list(tmp+'C')
	matrix[i] = tmp

matrix4f = np.delete(matrix4f,deleteIndex,0)
matrix0f = np.delete(matrix0f,deleteIndex,0)

matrix4f[0] = list(refSeq4f)
matrix0f[0] = list(refSeq0f)

df4f = pd.DataFrame(matrix4f)
# df4f = df4f.loc[:,df4f.iloc[0]!='N']   
# nunique = df4f.apply(pd.Series.nunique)
# cols_to_drop = nunique[nunique == 1].index
df4f=df4f.drop(cols_to_drop, axis=1)
df0f = pd.DataFrame(matrix4f)
df0f = df0f.loc[:,df0f.iloc[0]!='N']  

df4f = df4f.transpose()
df0f = df0f.transpose()

output=pd.DataFrame()
for index,r in df4f.iterrows():
	# print(j)
	div = 0
	af = 0
	
	ref = r.iloc[0]
	out = r.iloc[-1]
	
	if(out == 'N' or out == '-'):
		continue
	elif((out != ref)): 
		div = 1
		af = 0
	else:
		AA = out
		AN = r.shape[0]-2
		AC = r.iloc[1:(r.shape[0]-1)].value_counts()

		if(AA not in AC.index):
			continue
		else:
			AC = AC[AC!=AC[AA]]
			if(len(AC) == 0):
				af=0
			else:
				af=AC[0]/AN
	tmp = pd.DataFrame({'rawDerivedAllele':af,'div':div,'type':'4fold'},index=[id])
	tmp = tmp.reset_index(drop=True)
	output = pd.concat([output,tmp])
output=pd.DataFrame()
for index,r in df0f.iterrows():
	# print(j)
	div = 0
	af = 0
	
	ref = r.iloc[0]
	out = r.iloc[-1]
	
	if(out == 'N' or out == '-'):
		continue
	elif((out != ref)): 
		div = 1
		af = 0
	else:
		AA = out
		AN = r.shape[0]-2
		AC = r.iloc[1:(r.shape[0]-1)].value_counts()

		if(AA not in AC.index):
			continue
		else:
			AC = AC[AC!=AC[AA]]
			if(len(AC) == 0):
				af=0
			else:
				af=AC[0]/AN
	tmp = pd.DataFrame({'rawDerivedAllele':af,'div':div,'type':'0fold'},index=[id])
	tmp = tmp.reset_index(drop=True)
	output = pd.concat([output,tmp])

div = output.groupby(['type'])['div'].sum().reset_index()
# div = div.pivot_table(index=['type'],columns=['type'],values='div').reset_index()


daf = output[output['rawDerivedAllele']!=0]

bins = np.arange(0,1.05,0.05)
labels = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]

daf['categories'] = pd.cut(daf['rawDerivedAllele'],bins=bins,labels=labels)

sfs = daf.groupby(['type','categories']).count().reset_index()
sfs['rawDerivedAllele'] = sfs['rawDerivedAllele'].fillna(0).astype(int)
sfs = sfs.groupby(['type'])['rawDerivedAllele'].apply(list).reset_index()

# sfs['p0'] = sum(sfs[sfs['type']=='4fold']['rawDerivedAllele'].iloc[0])
# sfs['pi'] = sum(sfs[sfs['type']=='0fold']['rawDerivedAllele'].iloc[0])

print(time.time() - start)

# coordinates = '640959,644314,644437,645716,647989,648361,648479,648886,648991,649765,649934,650293,652786,653288,654247,654700,655274,655866,656712,657384,705871,707766,6973477,6973508,6973574,6973684,6973745,6975132,6975514,6975609,6975731,6975856,6975925,6976300,6976373,6981691,6981760,6983250,6983316,6983518,6983583,6985013,6987119,6987200,6999742,6999856'
# coordinates = array(coordinates.split(',')).astype(int).tolist()
# coordinates =  [coordinates[i:i+2] for i in range(0, len(coordinates), 2)]
