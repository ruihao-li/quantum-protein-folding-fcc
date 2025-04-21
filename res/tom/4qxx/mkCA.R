rm(list=ls())
setwd("/Users/raubenb/Desktop/gitrepos/protein-folding-qc/res/tom/4qxx")
system("grep CA pdbs/4qxx.pdb > pdbs/4qxxCA.pdb")
system("grep -v ANISOU pdbs/4qxxCA.pdb > pdbs/4qxx_CA.pdb")
system("rm pdbs/4qxxCA.pdb")# don't know how to use ANISOU lines, so remove file for now 
