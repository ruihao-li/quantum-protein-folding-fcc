rm(list=ls())
library(tidyverse) 
library(Rpdb)
(x <- atoms(recname = "ATOM", eleid = 1:2, elename = c("F"), alt = "",
           resname = "BCC", chainid = c("A"), resid = 1:2, insert = " ",
           x1 = c(0,1), x2 = c(0, 1), x3 = c(0,1), occ = c(0.0), temp = c(1.0),
           segid = ""))
#Guo ZY, Kraka E, Cremer D. Description of Local and Global Shape Properties of Protein Helices. J Mol Model. 2013;19(7):2901–2911.
# The radius, pitch and PPT of an α-helix are 2.3, 5.5 and 3.6 
# sc=2.3 
# 2.3*sqrt(3) # 3.98
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3892923/
# # distance between consecutive C α atoms are distributed normally with a mean of 3.8 Å and standard deviation of 0.04 Å.
# sc=2.194 # 3.8 A in AA space is sqrt(3)=1.7 in diamond space, so 3.8/sqrt(3)= 2.194 is scale factor
# 2.2*sqrt(3) # 3.810512
# # but Trp-cage specificaly yields 3.87 as the average distance between AA (removing first AA outlier)
# (sc=3.87/sqrt(3)) # 2.234346 is between 2.2 and 2.3, so use it or better yet
# sc=2.25 # keep things better rounded and not changing with each new protein (also 4x2.25=9 angstroms per unit cell edge)
# 3.8=sqrt(2)*sc
(sc=3.8/(sqrt(3))) # 2.193931
sc=2.2 # round to make files nicer
2.2*sqrt(3) #3.81
x$x1=sc*x$x1
x$x2=sc*x$x2
x$x3=sc*x$x3
x=pdb(x)
x$cryst1 <- cryst1(abc = c(2*sc, 2*sc,2*sc), abg = c(90,90,90), sgroup = "P1")
x$conect <- conect(x,safety=3.4) # a numeric value used to extend the atomic radii (=> radii are 0.2)
x$conect

visualize(x)
y <- replicate(x, a.ind = 0:12, b.ind = 0:12, c.ind = 0:12)
y$conect <- conect(y,safety=3.4) #maybe this works because we are extending radii of 0.25
visualize(y)
write.pdb(y,"trpLat/pdbs/BCC/LatBCC.pdb")
write_tsv(y$conect,"trpLat/pdbs/BCC/conBCC.txt",col_names=F)
