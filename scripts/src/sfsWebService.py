import time 
import textwrap
from pyfaidx import Fasta
start = time.time()


text = textwrap.dedent(tmp).strip()
file = Fasta('2l.seq')
coordinates = pd.read_csv('test.bed',sep='\t',header=None).iloc[:,[1,2]].values.tolist() 

samples = list(file.keys())

matrix = np.empty([len(samples),len(file.get_spliced_seq(samples[0], coordinates))+1],dtype='str')

# matrix[0] = file.get_spliced_seq(samples[0], coordinates).seq.upper()
degenCode = degenerate(file.get_spliced_seq(samples[0], coordinates).seq.upper()+'C')

deleteIndex = list()
seq = []

for i in range(2,len(samples),1):
	tmp = file.get_spliced_seq(samples[i],coordinates).seq.upper()

	if(tmp == ('N' * len(tmp))):
		deleteIndex.append(i)
	else:
		matrix[i] = list(tmp+'C')

matrix = np.delete(matrix,deleteIndex,0)
matrix[0] = list(degenCode);matrix[1] = list(file.get_spliced_seq(samples[0], coordinates).seq.upper()+'C')


start = time.time()
output = list()
for x in np.nditer(matrix,flags=['external_loop'], order='F'): 
	if((x[0] == '4') or (x[0] == '0')):
		
		if(x[0] == '4'):
			functionalClass = '4fold'
		else:
			functionalClass = '0fold'

		div = 0
		af = 0

		REF = x[1]
		AA = x[-1]

		if(out == 'N' or out == '-'):
			next
		elif((AA != REF) and (np.unique(x[2:-1]).shape[0] == 1)): 
			div = 1
			af = 0
		else:
			AN = x[2:-1].shape[0]
			AC = pd.DataFrame(data=np.unique(x[2:-1], return_counts=True)[1],index=np.unique(x[2:-1], return_counts=True)[0])
 			
			if(AA not in AC.index):
				next
			else:
				AC = AC[AC.index!=AA]
				if(len(AC) == 0):
					AF=0
				else:
					AF=AC.iloc[0]/AN
					AF = AF.iloc[0]
		# tmp = pd.DataFrame({'derivedAlleleFrequency':AF,'div':div,'type':functionalClass})
		# tmp = tmp.reset_index(drop=True)
		tmp = [AF,div,functionalClass]
		# print(tmp)
		output.append(tmp)
	else:
		next
time.time() - start
start = time.time()
output = pd.DataFrame()
for index,r in df.iterrows():
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
time.time()-start
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
