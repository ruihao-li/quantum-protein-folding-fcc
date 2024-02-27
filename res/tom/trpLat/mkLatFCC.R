rm(list=ls())
library(tidyverse) 
library(Rpdb)
(x <- atoms(recname = "ATOM", eleid = 1:4, elename = c("F"), alt = "",
           resname = "FCC", chainid = c("A"), resid = 1:4, insert = " ",
           x1 = c(0,1,0,1), x2 = c(0, 1, 1, 0), x3 = c(0,0,1,1), occ = c(0.0), temp = c(1.0),
           segid = ""))
(sc=3.8/(sqrt(2))) # 2.687006
x$x1=sc*x$x1
x$x2=sc*x$x2
x$x3=sc*x$x3
x=pdb(x)
x$cryst1 <- cryst1(abc = c(2*sc, 2*sc,2*sc), abg = c(90,90,90), sgroup = "P1")
2*2.687006 # 5.37 is length of unit cell
y <- replicate(x, a.ind = 0:10, b.ind = 0:10, c.ind = 0:10)
y$conect <- conect(y,safety=4)
visualize(y)
write.pdb(y,"trpLat/pdbs/FCC/LatFCC.pdb")
write_tsv(y$conect,"trpLat/pdbs/FCC/conFCC.txt",col_names=F)
