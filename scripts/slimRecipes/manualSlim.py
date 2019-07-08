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
args.rb=0
args.sb=0.1
args.sd=-0.2
args.dominanceCoef=0.5
args.burnin=10000
args.generations=21000
args.bins=20
args.output='alphaNeutral'
args.replica=1
args.path='/home/jmurga/mkt/201902/rawData/dmel/simulations/'


library(iMKT)
dafn <- fread('/home/jmurga/mkt/201902/rawData/dmel/simulations/wdFraction/alphaNeutral/alphaNeutraldaf0.tab'); divn <- fread('/home/jmurga/mkt/201902/rawData/dmel/simulations/wdFraction/alphaNeutral/alphaNeutraldiv0.tab')
n <- iMKT(daf=dafn[,c(1,2,3)],divergence=divn[,c(1,2,4,5)],xlow=0,xhigh=1,plot=F)


dafa <- fread('/home/jmurga/mkt/201902/rawData/dmel/simulations/wdFraction/adaptive1/adaptive1daf0.tab'); diva <- fread('/home/jmurga/mkt/201902/rawData/dmel/simulations/wdFraction/adaptive1/adaptive1div0.tab')
a <- iMKT(daf=dafa[,c(1,2,3)],divergence=diva[,c(1,2,4,5)],xlow=0,xhigh=1,plot=T)

daf <- fread('/home/jmurga/mkt/201902/rawData/dmel/simulations/wdFraction/adaptive2/adaptive2daf0.tab'); div <- fread('/home/jmurga/mkt/201902/rawData/dmel/simulations/wdFraction/adaptive2/adaptive2div0.tab')
a2 <- iMKT(daf=daf[,c(1,2,3)],divergence=div[,c(1,2,4,5)],xlow=0,xhigh=1,plot=T)

colnames(daf) <- c('daf','Pi','P0')
colnames(div) <- c('Di','D0','alpha')
div$mi <- 1e8
div$m0 <- 1e8

