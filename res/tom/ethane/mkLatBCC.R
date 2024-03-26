rm(list=ls())
library(tidyverse) 
library(Rpdb)
setwd("~/GH/protein-folding-qc/res/tom/ethane")
(x <- atoms(recname = "ATOM", eleid = 1:2, elename = c("F"), alt = "",
           resname = "BCC", chainid = c("A"), resid = 1:2, insert = " ",
           x1 = c(0,1), x2 = c(0, 1), x3 = c(0,1), occ = c(0.0), temp = c(1.0),
           segid = ""))

(sc=1.54/(sqrt(3))) 
x$x1=sc*x$x1
x$x2=sc*x$x2
x$x3=sc*x$x3
x=pdb(x)
x$cryst1 <- cryst1(abc = c(2*sc, 2*sc,2*sc), abg = c(90,90,90), sgroup = "P1")

# x$conect <- conect(x,safety=) # a numeric value used to extend the atomic radii (=> radii are 0.2)
# x$conect
# 
# visualize(x)
y <- replicate(x, a.ind = 0:1, b.ind = 0:1, c.ind = 0:1)
(y$conect <- conect(y,safety=1.5) )
# visualize(y)
write.pdb(y,"pdbs/LatBCC.pdb")
# write_tsv(y$conect,"pdbs/conBCC.txt",col_names=F)
