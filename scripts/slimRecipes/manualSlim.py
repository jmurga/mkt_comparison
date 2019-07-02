import argparse
import sys
import os
import random
import subprocess
import pandas as pd
# import rpy2
# import rpy2.robjects as robjects
from string import Template
from tempfile import NamedTemporaryFile


parser = argparse.ArgumentParser(description="SLiM simulations based on baseline recipe extract from asymptoticMK github, based on Messer et al 2013.")

args = parser.parse_args()

args.recipe='wdFraction'
args.length=1e7
args.mutRate=1e-9
args.recombRate=1e-7
args.ancSize=1000
args.rb=0.0005
args.sb=0
args.sd=0
args.dominanceCoef=0.5
args.burnin=10000
args.generations=210000
args.bins=20
args.output='alphaNeutral'
args.replica=10
args.path='/home/jmurga/mkt/201902/rawData/dmel/simulations/'