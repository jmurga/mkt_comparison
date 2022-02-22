#!/bin/bash

folder="/home/jmurga/mkt/201902/"


for f in theta_0.01 theta_0.0001 rho_0.01 rho_0.0001 alpha_0.1 alpha_0.7 base_weaky base_weakly_bgs
do
	mkdir -p ${folder}/raw_data/simulations/${f}/daf
	mkdir -p ${folder}/raw_data/simulations/${f}/div
done

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.001 -d theta=0.001 -d pposH=0.000213339 -d nGenes=7 -d nF={} -d p=\'${folder}/raw_data/simulations/base\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.001 -d theta=0.001 -d pposH=0.000213339 -d nGenes=1 -d nF={} -d p=\'${folder}/raw_data/simulations/genes_1k\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.001 -d theta=0.001 -d pposH=0.000213339 -d nGenes=14 -d nF={} -d p=\'${folder}/raw_data/simulations/genes_28k\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 



parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.001 -d theta=0.01 -d pposH=0.000213339 -d nGenes=7 -d nF={} -d p=\'${folder}/raw_data/simulations/theta_0.01\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.001 -d theta=0.0001 -d pposH=0.000213339 -d nGenes=7 -d nF={} -d p=\'${folder}/raw_data/simulations/theta_0.0001\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.01 -d theta=0.001 -d pposH=0.000213339 -d nGenes=7 -d nF={} -d p=\'${folder}/raw_data/simulations/rho_0.01\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.0001 -d theta=0.001 -d pposH=0.000213339 -d nGenes=7 -d nF={} -d p=\'${folder}/raw_data/simulations/rho_0.0001\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.001 -d theta=0.001 -d pposH=3.56199420564694e-05 -d nGenes=7 -d nF={} -d p=\'${folder}/raw_data/simulations/alpha_0.1\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 

parallel -j 20 -u "slim -d N=1000 -d negStrength=-2000 -d posStrength=250 -d shape=0.3 -d rho=0.001 -d theta=0.001 -d pposH=0.0007477524342020683 -d nGenes=7 -d nF={} -d p=\'${folder}/raw_data/simulations/alpha_0.7\' ${folder}/scripts/slim_recipes/mktComparison.slim" ::: `seq 1 2000` 


 parallel -j20 -u "slim -d sampleSize=20 -d codingLength=2000 -d muBgs=7.51e-19 -d N=1000 -d sS=500 -d wS=5 -d pposH=0.000093 -d pposL=0.006904 -d output=\'${folder}/raw_data/simulations/base_weakly/\' -d nF={} weakly.slim" ::: `seq 1 2000`

 parallel -j20 -u "slim -d sampleSize=20 -d codingLength=2000 -d muBgs=1.21e-06 -d N=1000 -d sS=500 -d wS=5 -d pposH=0.000093 -d pposL=0.006904 -d output=\'${folder}/raw_data/simulations/base_weakly_bgs/\' -d nF={} weakly.slim" ::: `seq 1 2000`
