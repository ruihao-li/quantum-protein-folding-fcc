rm(list=ls())
library(tidyverse) 
library(Rpdb)
(x <- atoms(recname = "ATOM", eleid = 1, elename = c("F"), alt = "",
           resname = "CUB", chainid = c("A"), resid = 1, insert = " ",
           # x1 = c(0.5), x2 = c(0.5), x3 = c(0.5), occ = c(0.0), temp = c(1.0),
           x1 = c(0), x2 = c(0), x3 = c(0), occ = c(0.0), temp = c(1.0),
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
2.25*sqrt(3)# 3.9
(sc=3.8) # distance between Calphas
3.8*sqrt(3) # 3-barrel pitch of 6.6 is closer to 5.5 than 9.0
x$x1=sc*x$x1
x$x2=sc*x$x2
x$x3=sc*x$x3
x=pdb(x)
# x$conect <- conect(x,safety=4)
# x$conect  
x$cryst1 <- cryst1(abc = c(sc, sc,sc), abg = c(90,90,90), sgroup = "P1")
# write.pdb(x,"trpDia/pdbs/LatCu1.pdb")
# before 9A x 6 was 54 A across 6x6x6 cube, here go with 15*3.8=57 at end
# y <- replicate(x, a.ind = 0:6, b.ind = 0:6, c.ind = 0:6)
y <- replicate(x, a.ind = 0:15, b.ind = 0:15, c.ind = 0:15)
y$conect <- conect(y,safety=4)
visualize(y)
write.pdb(y,"pdbs/CUB/LatCUB.pdb")
write_tsv(y$conect,"pdbs/CUB/conCUB.txt",col_names=F)
