# Comparison of four hierachical McDonald and Kreitman test approaches using real and simulated population data

Here, we perform a comparison of four hierarchical MKT methods: (i) the standard (original) MKT (sMKT, McDonald and Kreitman 1991); (ii) the Fay, Wickoff and Wu correction (FWWMKT, Fay 2001); (iii) a new imputation method based on FWWMKT (imputationMKT)the Extended MKT (eMKT, Mackay et al. 2012); and (iv) the asymptotic MKT (aMKT, Messer and Petrov 2012). 

Both real and simulated population genomics data are used to assess their performance for different evolutionary scenarios. Genome-wide DNA variation data come from two *Drosophila melanogaster* and two human populations (one colonizing and one ancestral population, in each case). Simulated data was generated with the SLiM 3 evolutionary framework (Haller ‎2019). We test several conditions including gene-to-gene vs gene concatenating analysis, or the effect of recombination on the power and bias of selection estimates of the different MKT methods.


This repository only include raw code to get main results. notebooks/ folder include Jupyter Notebooks running on Python 3.6 kernel to execute step by step the pipeline. src/ folder contain raw scripts to needed to execute the pipelin. Please note that multiple step could be parallelized, in this case create yourself customs bash scripts or run it on your server manually.  

Pipeline were developed with the following [conda enviroment](https://github.com/jmurga/mkt_comparison/blob/master/imp.yml) in local server: 100GB RAM and 16 Intel(R) Xeon(R) CPU.  

Pipelines execution requiere to download  files. Paths would need to be changed too.
