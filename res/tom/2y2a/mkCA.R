rm(list=ls())
setwd("/Users/raubenb/Desktop/gitrepos/protein-folding-qc/res/tom/2y2a")
system("grep CA pdbs/2y2a.pdb > pdbs/2y2aCA.pdb")
system("grep -v ANISOU pdbs/2y2aCA.pdb > pdbs/2y2a_CA.pdb")
system("rm pdbs/2y2aCA.pdb")# don't know how to use ANISOU lines, so remove file for now 
