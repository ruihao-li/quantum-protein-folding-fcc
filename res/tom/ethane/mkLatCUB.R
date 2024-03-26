rm(list=ls())
library(tidyverse) 
library(Rpdb)
setwd("~/GH/protein-folding-qc/res/tom/ethane")
(x <- atoms(recname = "ATOM", eleid = 1, elename = c("F"), alt = "",
           resname = "CUB", chainid = c("A"), resid = 1, insert = " ",
           x1 = c(0), x2 = c(0), x3 = c(0), occ = c(0.0), temp = c(1.0),
           segid = ""))
(sc=1.54) # C-C bond length 
x$x1=sc*x$x1
x$x2=sc*x$x2
x$x3=sc*x$x3
x=pdb(x)
x$cryst1 <- cryst1(abc = c(sc, sc,sc), abg = c(90,90,90), sgroup = "P1")
y <- replicate(x, a.ind = 0:2, b.ind = 0:2, c.ind = 0:2)
y$conect <- conect(y,safety=1.6)
# visualize(y)
write.pdb(y,"pdbs/LatCUB.pdb")
write_tsv(y$conect,"pdbs/conCUB.txt",col_names=F)
