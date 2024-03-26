rm(list=ls())
library(tidyverse) 
library(Rpdb)
setwd("~/GH/protein-folding-qc/res/tom/ethane")
x <- atoms(recname = "ATOM", eleid = 1:8, elename = c("N","F"), alt = "",
           resname = "DIA", chainid = c("A","B"), resid = 1:8, insert = " ",
           x1 = c(0,1,0,1,2,3,2,3), x2 = c(0,1,2,3,0,1,2,3), x3 = c(0,1,2,3,2,3,0,1), occ = c(0.0,0.0), temp = c(1.0,1.0),
           segid = "")
sc=1.54/sqrt(3)
x$x1=sc*x$x1
x$x2=sc*x$x2
x$x3=sc*x$x3
x=pdb(x)
x$conect <- conect(x,safety=1.5)
x$conect  
x$cryst1 <- cryst1(abc = c(4*sc, 4*sc, 4*sc), abg = c(90,90,90), sgroup = "P1")

y <- replicate(x, a.ind = 0:1, b.ind = 0:1, c.ind = 0:1)
y$conect <- conect(y,safety=1.5)
visualize(y)

write.pdb(y,"pdbs/LatDIA.pdb")
write_tsv(y$conect,"pdbs/conDIA.txt",col_names=F)
