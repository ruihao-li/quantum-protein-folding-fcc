rm(list=ls())
setwd("~/GH/protein-folding-qc/res/tom/chigLat")
system("grep CA pdbs/5awl.pdb > pdbs/5awlCA.pdb")
system("grep -v ANISOU pdbs/5awlCA.pdb > pdbs/5awl_CA.pdb")
system("rm pdbs/5awlCA.pdb")# don't know how to use ANISOU lines, so remove file for now 
