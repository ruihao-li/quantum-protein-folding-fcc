rm(list=ls())
library(tidyverse) 
library(Rpdb)
setwd("~/GH/protein-folding-qc/res/tom/ethane")
(x <- atoms(recname = "ATOM", eleid = 1:4, elename = c("F"), alt = "",
           resname = "FCC", chainid = c("A"), resid = 1:4, insert = " ",
           x1 = c(0,1,0,1), x2 = c(0, 1, 1, 0), x3 = c(0,0,1,1), occ = c(0.0), temp = c(1.0),
           segid = ""))
(sc=1.54/(sqrt(2)))
x$x1=sc*x$x1
x$x2=sc*x$x2
x$x3=sc*x$x3
x=pdb(x)
x$cryst1 <- cryst1(abc = c(2*sc, 2*sc,2*sc), abg = c(90,90,90), sgroup = "P1")
y <- replicate(x, a.ind = 0:1, b.ind = 0:1, c.ind = 0:1)
(y$conect <- conect(y,safety=1.36))
# visualize(y)
write.pdb(y,"pdbs/LatFCC.pdb")
# write_tsv(y$conect,"pdbs/conFCC.txt",col_names=F)
