#!/home/jmurga/.conda/envs/simulations/bin/python
######################################################################
###############SCRIPT TO AUTOMATIZE SLIM SIMULATIONS##################
#####################CREATES TEMPORAL RECIPES#########################
##############BASED ON THE SELECTED RECIPE AND ARGPARSE MODULE########
########MODIFIED FROM https://github.com/bodkan/nea-over-time#########
######################################################################
import argparse
import sys
import os
import random
import subprocess
import pandas as pd
import numpy as np
# import rpy2
# import rpy2.robjects as robjects
from string import Template
from tempfile import NamedTemporaryFile
# from rpy2 import robjects
# from rpy2.robjects import pandas2ri

# pandas2ri.activate()


# Write lists to an .deleteriousFreqata file. Two functions in odeleteriousFreqer to save differents global variable name.
def savedeleteriousFreqata(ls,filename,path,variableName):
	robjects.r.assign('{}', robjects.Vector(ls)).format(variableName)
	robjects.r("save({}, file='{}{}')".format(variableName,path,filename))


def slimVector(xs):
	"""Convert a list of numbers to a SLiM code creating the same list."""
	return "c(" + ",".join(str(x) for x in xs) + ")"

if __name__ == "__main__":
	'''Parse arguments and show the required inputs if only name is given to command line'''
	parser = argparse.ArgumentParser(description="SLiM simulations based on baseline recipe extract from asymptoticMK github, based on Messer et al 2013.")

	# Required arguments
	parser.add_argument("--recipe", type = str, required = True, help = "Slim recipe to execute")
	parser.add_argument("--length", type = float, required = 1e7, help = "0-based start and stop SLiM length to simulate.")
	parser.add_argument("--mutRate", type = float, required = 1e-9, help = "Mutation rate in the simulated region")
	parser.add_argument("--recombRate", type = float, required = 1e-7, help = "Recombination rate in the simulated region")
	parser.add_argument("--dominanceCoef", type = float, required = 0.5, help = "")
	parser.add_argument("--neutralFreq", type = float, required = 0.5, help = "")
	parser.add_argument("--deleteriousFreq", type = float, required = 0.5, help = "")
	parser.add_argument("--beneficialFreq", type = float, required = 0.0005, help = "")
	parser.add_argument("--deleteriousFitness", type = float, required = 0.1, help = "")
	parser.add_argument("--beneficialFitness", type = float, required = 0.2, help = "")
	parser.add_argument("--ancSize", type = int, required = 1000, help = "Effective population size of the ancestral population")
	parser.add_argument("--burnin", type = int, required = 10000, help = "Burnin period")
	parser.add_argument("--generations", type = int, required = 210000, help = "Burnin period")
	parser.add_argument("--bins", type = int, required = 20, help = "Burnin period")
	parser.add_argument("--replica", type = int, required = 20, help = "Number of replica by scenario")

	parser.add_argument("--scenario", type = str,help = "Output name without extension")
	parser.add_argument("--output", type = str,help = "output")
	
	# Default arguments
	parser.add_argument("--path", type = str, default = '/home/jmurga/mkt/201902/rawData/simulations', help = "Path to output file")

	
	# Parsing common arguments
	args = parser.parse_args()

	output = args.path + '/' + args.recipe + '/' + args.scenario + '/' + args.output
	regions = pd.DataFrame({'start':0, 'end':args.length}, index=[0])
	regions = regions[['start','end']]
	burnin = 10 * args.ancSize

	# Selecting baseline recipe

	slimRecipe = Template(open("/home/jmurga/mkt/201902/scripts/slimRecipes/" + args.recipe + '.slim', "r").read())

	mapping = {
		# 'genomicElements' : genomicElements,
		'mutRate' : args.mutRate,
		'length' : int(args.length),
		'recombRate' : args.recombRate,
		'neutralFreq' : float(args.neutralFreq),
		'beneficialFreq' : float(args.beneficialFreq),
		'deleteriousFreq' : float(args.deleteriousFreq),
		'beneficialFitness' : float(args.beneficialFitness),
		'deleteriousFitness' : float(args.deleteriousFitness),
		'h' : float(args.dominanceCoef),
		'ancSize' : args.ancSize,
		'generations' : args.generations,
		'burnin' : args.burnin,
		'bins' : int(args.bins),
		'output' : output
	}



	# print(slimRecipe.substitute(mapping))
	
	with NamedTemporaryFile("w") as slim_file:
		print(slimRecipe.substitute(mapping), file= slim_file,flush=True)

		# Create directory
		os.makedirs(output,exist_ok=True)

		for i in range(0,args.replica,1):
			print(i)

			# Opening slim procces and save custom string output in python variable v
			# str(random.randint(1, 10**13))
			slimResults = subprocess.run(["/home/jmurga/mkt/201902/software/SLiM3.3/slim", "-s", str(random.randint(1, 10**13)),slim_file.name],universal_newlines=True,stdout=subprocess.PIPE)
			
			# Parsing string output, we checked position on slim custom printed output and procces each variable taking into account correspondent positions. Excluding recipe execution info
			slimResults = slimResults.stdout.split('\n')
			
			# Need to extract position [0] on list due to split function return a list based on split pattern. Each line is an element at the list
			if(args.recipe=='baseline'):

				rawD0 = float(slimResults[15:16][0].split(':')[1])
				rawD = float(slimResults[16:17][0].split(':')[1])
				rawTrueAlpha = float(slimResults[17:18][0].split(':')[1])
	
				rawDaf = slimResults[18:]

				rawDaf = [x.split('\t') for x in rawDaf] 
				rawDaf = rawDaf[:-1]
				header = rawDaf[0]
				rawDaf = rawDaf[1:]

				daf = pd.DataFrame(rawDaf,columns=header)

				# Save results in lists based on scenario replicas
				# d0.append(rawD0)
				# d.append(rawD)
				# trueAlpha.append(rawTrueAlpha)
				# divergenceAndTrueAlpha.append(div)
				# listDaf.append(daf)

				print(daf)

				div = pd.DataFrame({'d':rawD,'d0':rawD0,'trueAlpha':rawTrueAlpha},index=[0])
				
				div.to_csv(args.path + args.recipe + '/' + args.output + '/' + args.output + 'div' + str(i)+'.tab',index=False,header=True, sep='\t')

				daf.to_csv(args.path + args.recipe + '/' + args.output + '/' + args.output+'daf'+str(i)+'.tab',index=False,header=True, sep='\t')

			elif(args.recipe=='wdEstimation'):
					
				# Cleaning results
				# Daf
				slimDaf = slimResults[slimResults.index('daf\tPi\tP0\tPneu\tPwd'):slimResults.index('D0\tDi\tm0\tmi\ttrueAlpha\tb')]
				slimDaf = [x.split('\t') for x in slimDaf]

				# Header
				h = slimDaf[0]

				# Daf
				d = slimDaf[1:]

				daf = pd.DataFrame(d,columns=h)
				print(daf)
				
				# Processing div info
				slimDiv = slimResults[slimResults.index('D0\tDi\tm0\tmi\ttrueAlpha\tb'):-1]
				slimDiv = [x.split('\t') for x in slimDiv]
				
				# Header
				h = slimDiv[0]
				# Divergence
				d = slimDiv[1]

				div = pd.DataFrame([d],columns=h)

				print(div)

				daf.to_csv(output + '/' + 'daf' + str(i) + '.tab',index=False,header=True, sep='\t')
				div.to_csv(output + '/' + 'div' + str(i)+'.tab',index=False,header=True, sep='\t')

		# outputDiv = args.path + args.output + '/'
		# outputDaf = args.path + args.output + '/'

		# Save as deleteriousFreqata object including all the scenarios
		# savedeleteriousFreqata(divergenceAndTrueAlpha,'div.deleteriousFreqata',outputDiv,div)
		# savedeleteriousFreqata(listDaf,'daf.deleteriousFreqata',outputDaf,daf)
