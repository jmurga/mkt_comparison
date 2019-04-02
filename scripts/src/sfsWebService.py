import time
import numpy as np
from pyfaidx import Fasta
def degenerancy(data):
	
	#DEGENERANCY DICTIONARIES
	degenerateCodonTable = {"TTC":'002',"TTT":'002',#PHE 002
						"TTA":'002',"TTG":'002',#LEU 002
						"CTT":'204',"CTC":'204',"CTA":'204',"CTG":'204',#LUE 204
						"ATG":'000',#MET 000
						"ATT":'003',"ATC":'003',"ATA":'003',#ILE 003
						"GTC":'004',"GTT":'004',"GTA":'004',"GTG":'004',#VAL 004
						"TCT":'004',"TCA":'004',"TCC":'004',"TCG":'004',#SER 004
						"CCT":'004',"CCA":'004',"CCC":'004',"CCG":'004',#PRO 004
						"ACT":'004',"ACA":'004',"ACC":'004',"ACG":'004',#THR 004
						"GCT":'004',"GCA":'004',"GCC":'004',"GCG":'004',#ALA 004
						"TAT":'002',"TAC":'002',#TYR 002
						"TAA":'000',"TAG":'000',"TGA":'000',#STOP 000
						"CAT":'002',"CAC":'002',#HIS 002
						"CAA":'002',"CAG":'002',#GLN 002
						"AAT":'002',"AAC":'002',#ASN 002
						"AAG":'002',"AAA":'002',#LYS 002
						"GAT":'002',"GAC":'002',#ASP 002
						"GAA":'002',"GAG":'002',#GLU 002
						"TGT":'002',"TGC":'002',#CYS 002
						"TGG":'000',#TRP 000
						"CGT":'204',"CGG":'204',"CGA":'204',"CGC":'204',#ARG 204
						"AGA":'002',"AGG":'002',#ARG 002
						"AGT":'002',"AGC":'002',#SER 002
						"GGT":'004',"GGA":'004',"GGC":'004',"GGG":'004',#GLY 004					
	}
	
	degenerancy = ''
	for i in range(0, len(data),3):
		codon = data[i:i+3]
		if(codon == 'NNN'):
			degenerancy += codon
		else:
			degenerancy += degenerateCodonTable[codon]

	return(degenerancy)
start = time.time()

file = Fasta('kilobase.fa')
# file = Fasta('kilobase.fa')
samples = list(file.keys())

matrix = np.empty([len(samples),len(file[samples[0]][:-1].seq)],dtype='str')
print(time.time() - start)

# Degenerate reference sequences (always first entry at file)
degenCode = degenerancy(file[samples[0]][:-1].seq.upper())

# List to append indexes if whole sequence at any population is len(seq) * 'N'
deleteIndex = list()
for i in range(1,len(samples),1):
	# Extract each sample sequence
	tmp = file[samples[i]][:-1].seq.upper()
	if(tmp == ('N' * len(tmp))):
		deleteIndex.append(i)
	else:
		matrix[i] = list(tmp)
print(time.time() - start)

# Delete lines
matrix = np.delete(matrix,deleteIndex,0)
# Put degenerancy in first ndarray element
matrix[0] = list(degenCode)
# NEED TO SOLVE THIS
solve = time.time()
d = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')  
print(time.time() - solve)


output = list()
print(time.time() - start)

for x in np.nditer(d, order='F',flags=['external_loop']): 
	print(x)
	REF = x[1]
	AA = x[-1]
	# Undefined Ancestra Allele. Try to clean out of the loop
	if(AA == 'N' or AA == '-'):
		next
	# Monomorphic sites. Try to clean out of the loop
	elif(np.unique(x[1:][np.where(x[1:]!='N')]).shape[0] == 1):
		next
	else:
		# print(np.unique(x[1:][np.where(x[1:]!='N')]).shape[0])
		pol = x[2:-1][np.where(x[2:-1]!='N')]
		if(x[0] == '4'):
			functionalClass = '4fold'
		else:
			functionalClass = '0fold'
		# CHECK IF POL != AA
		if((AA != REF) and (np.unique(pol).shape[0] == 1) and (np.unique(pol)[0] != AA)):
		# if((AA != REF)): 
			div = 1; AF = 0
			tmp = [AF,div,functionalClass]
			output.append(tmp)
		else:
			AN = x[2:-1].shape[0]
			AC = pd.DataFrame(data=np.unique(x[2:-1], return_counts=True)[1],index=np.unique(x[2:-1], return_counts=True)[0])
			div = 0
			if(AA not in AC.index):
				next
			else:
				AC = AC[AC.index!=AA]
				if(len(AC) == 0):
					next
				else:
					AF=AC.iloc[0]/AN
					AF = AF.iloc[0]
			tmp = [AF,div,functionalClass]
			output.append(tmp)
print(time.time() - start)

df = pd.DataFrame(output)
df['id'] = 'uploaded'
df.columns = ['derivedAlleleFrequency','d','functionalClass','id']

# Extract divergence data
div = df[['derivedAlleleFrequency','id','functionalClass','d']]
div = div.groupby(['id','functionalClass'])['d'].count().reset_index()
div = div.pivot_table(index=['id'],columns=['functionalClass'],values='d').reset_index()
div = div[['0fold','4fold']]
div['mi'] =  matrix[0][np.where(matrix[0]=='0')].shape
div['m0'] =  matrix[0][np.where(matrix[0]=='4')].shape     
div.columns = ['Di','D0','mi','m0']
# div = div.pivot_table(index=['functionalClass'],columns=['functionalClass'],values='div').reset_index()

# Create SFS pd.DataFrame by functionClass and 20 frequency bin
daf = df[df['d']!=1][['derivedAlleleFrequency','functionalClass','id']]
bins = np.arange(0,1.05,0.05)
labels = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]
daf['categories'] = pd.cut(daf['derivedAlleleFrequency'],bins=bins,labels=labels)
daf = daf.groupby(['functionalClass','id','categories']).count().reset_index()
sfs = pd.DataFrame({'daf':daf['categories'].unique(),'P0':daf[daf['functionalClass']=='4fold']['derivedAlleleFrequency'].reset_index(drop=True),'Pi':daf[daf['functionalClass']=='0fold']['derivedAlleleFrequency'].reset_index(drop=True)})

print(time.time() - start)
