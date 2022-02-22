import numpy as np
import pandas as pd
import pyfaidx as px

def reverseComplement(seq):
	# Sequence on dgn/sp contains M. No idea why
	seq = seq[::-1]
	letters = list(seq)
	baseComplement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
#	letters = [baseComplement[base] for base in letters]
	letters = [baseComplement[base] if base in baseComplement.keys() else base for base in letters]
	seq = ''.join(letters)

	return(seq)


def nDistribution(data,population,FASTAS='/data/shared/dgn'):
	''' Function to retrieve the number of undetermined nucleotide by largest transcript positions. Need as input a variable with transcript and coordinates and the population name to extract from dgn folder. '''
	for index, row in data.iterrows():    
		print(row['id'])
		# Convert CDS list into numeric array
		coordinates = np.array(row['coordinates'].split(',')).astype(int).tolist()
		coordinates =  [coordinates[i:i+2] for i in range(0, len(coordinates), 2)]
		# Open ref and outgroup
		ref = px.Fasta(FASTAS + '/ref/Chr' + row['chr'] +'.fasta')  
		## Extract ref and outgroup seq
		refSeq = ref.get_spliced_seq(row['chr'],coordinates).seq.upper()    
		if('N' not in refSeq):
			# Open population multifasta
			popFasta = px.Fasta(FASTAS + '/alignments/' + population + '_Chr' + row['chr'] +'.seq')
			#Extract samples
			samples = list(popFasta.keys())
			matrixDna = np.empty([len(samples)+1,len(refSeq)],dtype='str')
			if(row['strand'] == '-'):            
				refSeq = reverseComplement(refSeq)
				matrixDna[0] = list(refSeq)
				for i in range(0,len(samples),1):
					tmp = popFasta.get_spliced_seq(samples[i], coordinates).seq.upper()
					tmp = reverseComplement(tmp)
					matrixDna[i+1] = list(tmp)
			else:
				matrixDna[0] = list(refSeq)
				for i in range(0,len(samples),1):
					tmp = popFasta.get_spliced_seq(samples[i], coordinates).seq.upper()
					matrixDna[i+1] = list(tmp)
			# Count occurences
			df = pd.DataFrame(matrixDna).transpose()
			for i,r in df.iterrows():
				if('N' in r.values):
					tmp = pd.DataFrame({'m':r.value_counts()['N']},index=[0])
					tmp.to_csv('/home/jmurga/mktComparison/rawData/dmel/ncall'+population+'.tab',header=False,index=False,mode='a')
				else:
					continue
		else:
			continue