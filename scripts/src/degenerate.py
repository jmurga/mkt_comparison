#!/usr/bin/python
from __future__ import division
from Bio import SeqIO
import os
import numpy as np
import pandas as pd
import sys
import heapq
import uuid
import copy


def reverse(dna):
  return dna[::-1]

def complement (dna):
   basecomplement={'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a','n':'n','-':'-'}
   letters=list(dna)
   letters=[basecomplement[base] for base in letters]
   return "".join(letters)

def reverse_complement(dna):
   dna=reverse(dna)
   dna=complement(dna)
   return dna

arg1 = sys.argv[1]#TFBS
arg2= sys.argv[2]#4fold file
arg3 = sys.argv[3]#0fold file

###############################################################
##################FUNCTIONS####################################
###############################################################
###############################################################
##################PARSING FILE WITH SEQ.IO#####################
###############################################################

def transpose(arg1):
	try:
		f = open(arg1)
	except IOError:                     
		print("The file, %s, does not exist" % arg1)
		return

	ids=[]
	allSeqs=[]
	seq_daf=[]
	seq_maf=[]

	# print(SeqIO.parse(arg1, """fasta"""))
	print("Parsing fasta file")
	for seq_record in SeqIO.parse(arg1, """fasta"""):

		if (seq_record.seq[0:0+3]=='TTA') or (seq_record.seq[0:0+3]=='TCA') or (seq_record.seq[0:0+3]=='CTA'):
			seq_record.seq=reverse_complement(seq_record.seq)
			allSeqs.append(seq_record.seq) #.seq=sequences, can save in other variable the IDs with .id at the iteration
			ids.append(seq_record.id)
			# print('>'+seq_record.id)
			# print(seq_record.seq)
		elif (seq_record.seq[0:0+3]=='ATG'):
			allSeqs.append(seq_record.seq) #.seq=sequences, can save in other variable the IDs with .id at the iteration
			ids.append(seq_record.id)
			# print(seq_record.seq)
		else:
			continue
	# print('+++++++in matrix+++++++++')
	print("Creating sequence matrix")

	matrix_sequences=np.empty([len(allSeqs),len(allSeqs[0])],dtype='str')

	for i in range(0,len(allSeqs),1):
		matrix_sequences[i]=list(allSeqs[i]) #transform 
	return(matrix_sequences,ids)



###############################################################
##################DEGENERATE SEQUENCES#########################
###############################################################
def degenerate(arg,argid,file4fold,file0fold):
	codon_table_4f={
		"TTC":'NNN',"TTT":'NNN',#PHE 002
	    "TTA":'NNN',"TTG":'NNN',#LEU 202
	   	"CTT":'NNT',"CTC":'NNC',"CTA":'NNA',"CTG":'NNG',#LUE 204
	   	"ATG":'NNN',#MET 000
	   	"ATT":'NNN',"ATC":'NNN',"ATA":'NNN',#ILE 003
		"GTC":'NNC',"GTT":'NNT',"GTA":'NNA',"GTG":'NNG',#VAL 004
		"TCT":'NNT',"TCA":'NNA',"TCC":'NNC',"TCG":'NNG',#SER 004
		"CCT":'NNT',"CCA":'NNA',"CCC":'NNC',"CCG":'NNG',#PRO 004
		"ACT":'NNT',"ACA":'NNA',"ACC":'NNC',"ACG":'NNG',#THR 004   
		"GCT":'NNT',"GCA":'NNA',"GCC":'NNC',"GCG":'NNG',#ALA 004
		"TAT":'NNN',"TAC":'NNN',#TYR 002
		"TAA":'NNN',"TAG":'NNN',"TGA":'NNN',#STOP 000
		"CAT":'NNN',"CAC":'NNN',#HIS 002
		"CAA":'NNN',"CAG":'NNN',#GLN 002
		"AAT":'NNN',"AAC":'NNN',#ASN 002
		"AAG":'NNN',"AAA":'NNN',#LYS 002
		"GAT":'NNN',"GAC":'NNN',#ASP 002
		"GAA":'NNN',"GAG":'NNN',#GLU 002
		"TGT":'NNN',"TGC":'NNN',#CYS 002 
		"TGG":'NNN',#TRP 000
		"CGT":'NNT',"CGG":'NNG',"CGA":'NNA',"CGC":'NNC',#ARG 204
		"AGA":'NNN',"AGG":'NNN',#ARG 002
		"AGT":'NNN',"AGC":'NNN',#SER 002
		"GGT":'NNT',"GGA":'NNA',"GGC":'NNC',"GGG":'NNG',#GLY 004
		"NNN":'NNN',
	   	"---":'---',
	   }
	codon_table_0f={
	  	"TTC":'TTN',"TTT":'TTN',#PHE 002
		"TTA":'NTN',"TTG":'NTN',#LEU 202
		"CTT":'NTN',"CTC":'NTN',"CTA":'NTN',"CTG":'NTN',#LUE 204
		"ATG":'ATG',#MET 000
		"ATT":'ATN',"ATC":'ATN',"ATA":'ATN',#ILE 003
		"GTC":'GTN',"GTT":'GTN',"GTA":'GTN',"GTG":'GTN',#VAL 004
		"TCT":'TCN',"TCA":'TCN',"TCC":'TCN',"TCG":'TCN',#SER 004
		"CCT":'CCN',"CCA":'CCN',"CCC":'CCN',"CCG":'CCN',#PRO 004
		"ACT":'ACN',"ACA":'ACN',"ACC":'ACN',"ACG":'ACN',#THR 004
		"GCT":'GCN',"GCA":'GCN',"GCC":'GCN',"GCG":'GCN',#ALA 004
		"TAT":'TAN',"TAC":'TAN',#TYR 002
		"TAA":'NNN',"TAG":'NNN',"TGA":'NNN',#STOP 000
		"CAT":'CAN',"CAC":'CAN',#HIS 002
		"CAA":'CAN',"CAG":'CAN',#GLN 002
		"AAT":'AAN',"AAC":'AAN',#ASN 002
		"AAG":'AAN',"AAA":'AAN',#LYS 002
		"GAT":'GAN',"GAC":'GAN',#ASP 002
		"GAA":'GAN',"GAG":'GAN',#GLU 002
		"TGT":'TGN',"TGC":'TGN',#CYS 002
		"TGG":'TGG',#TRP 000
		"CGT":'NGN',"CGG":'NGN',"CGA":'NGN',"CGC":'NGN',#ARG 204
		"AGA":'AGN',"AGG":'AGN',#ARG 002
		"AGT":'AGN',"AGC":'AGN',#SER 002
		"GGT":'GGN',"GGA":'GGN',"GGC":'GGN',"GGG":'GGN',#GLY 004
		"NNN":'NNN',
               }

	matrix_4fold=matrix_sequences
	matrix_0fold=copy.deepcopy(matrix_sequences)

	degenerate_0fold=''
	degenerate_4fold=''
	# print(matrix_4fold)
	# rev_0fold=''
	# rev_4fold=''

	myfile4fold = open(file4fold, 'w')
	myfile0fold = open(file0fold, 'w')

	dna=arg[0].tostring()
 	# print(len(dna))
 	# print(dna)
	for i in range(0, len(dna),3):
		codon=dna[i:i+3]
		# print(codon)
		if(dna[0:0+3]=='ATG') or (dna[0:0+3]=='ATN') or (dna[0:0+3]=='ANG') or (dna[0:0+3]=='NTG') or (dna[0:0+3]=='ATA') or (dna[0:0+3]=='---'):
			if ("N" not in codon) and ("-" not in codon):
				degenerate_4fold+=codon_table_4f[codon]
				degenerate_0fold+=codon_table_0f[codon]
			else:
				codon='NNN'
				degenerate_4fold+=codon
				degenerate_0fold+=codon

		# elif (dna[0:0+3]=='TTA') or (dna[0:0+3]=='TCA') or (dna[0:0+3]=='CTA'):
		# 	negative=reverse_complement(dna)
		# 	codon_negative=negative[i:i+3]
		# 	if ("N" not in codon_negative) and ("-" not in codon_negative):
		# 			rev_4fold+=codon_table_4f[codon_negative]
		# 			rev_0fold+=codon_table_0f[codon_negative]
  # print(rev)
  	# print(degenerate_4fold,degenerate_0fold)
  	# print(matrix_4fold.size,matrix_sequences.size)
	# if (rev_4fold == '') and (rev_0fold == ''):
	# print(arg.shape,arg.shape[1])

	for j in range(0,arg.shape[1],1):
		# print(j)
		if(degenerate_4fold[j]=='N'):
			matrix_4fold[:,j]='N'  
		if(degenerate_0fold[j]=='N'):
			matrix_0fold[:,j]='N'
	# print(matrix_4fold)

	for row in range(0,arg.shape[0],1):
		myfile4fold.write('>'+argid[row]+'\n')
		myfile4fold.write(matrix_4fold[row].tostring()+'\n')

		myfile0fold.write('>'+argid[row]+'\n')
		myfile0fold.write(matrix_0fold[row].tostring()+'\n')

	# else:
	# 	for j in range(0,arg.shape[1],1):
	# 		if(rev_4fold[j]=='N'):
	# 			matrix_4fold[:,j]='N'  
	# 		elif (rev_0fold[j]=='N'):
	# 			matrix_0fold[:,j]='N'

	# 	for row in range(0,arg.shape[0],1):
	# 		myfile4fold.write('>'+argid[row]+'\n')
	# 		myfile4fold.write(matrix_4fold[row].tostring()+'\n')

	# 		myfile0fold.write('>'+argid[row]+'\n')
	# # 		myfile0fold.write(matrix_0fold[row].tostring()+'\n')

	# myfile4fold.close()
	# myfile0fold.close()

###############################################################
##################COMMAND LINES################################
###############################################################
print('++++++++++Creating 4fold and 0fold file++++++++++')
matrix_sequences,ids=transpose(arg1)
# print(matrix_sequences,ids)
degenerate(matrix_sequences,ids,arg2,arg3)
# 

