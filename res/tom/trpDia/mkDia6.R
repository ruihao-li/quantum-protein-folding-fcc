rm(list=ls())
library(tidyverse) # for things like mutate()
library(bio3d) # to read and write pdbs
pdb <- read.pdb("trpDia/pdbs/diamond.pdb") # cubic unit cell four units across
(p0=pdb$atom)               # pdb in angstrom. Start with 4 angstroms across repeated cell
h=1
N=6
L=vector(mode="list",length=N^3)
for (i in 0:(N-1))
  for (j in 0:(N-1))
    for (k in 0:(N-1)) {
      L[[h]]=p0%>%mutate(eleno=eleno+8*(h-1),resno=resno+8*(h-1),x=x+4*i,y=y+4*j,z=z+4*k)
      h=h+1  
    }
tail(p<-bind_rows(L))

myf=function(d) matrix(as.vector(t(as.matrix(d%>%select(x:z)))),nrow=1)
nL=lapply(L,myf)
xyz=matrix(unlist(nL),nrow=1)
pdb$atom=p
pdb$xyz=xyz
#Guo ZY, Kraka E, Cremer D. Description of Local and Global Shape Properties of Protein Helices. J Mol Model. 2013;19(7):2901–2911.
# The radius, pitch and PPT of an α-helix are 2.3, 5.5 and 3.6 
sc=2.3 
2.3*sqrt(3) # 3.98
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3892923/
# distance between consecutive C α atoms are distributed normally with a mean of 3.8 Å and standard deviation of 0.04 Å.
sc=2.194 # 3.8 A in AA space is sqrt(3)=1.7 in diamond space, so 3.8/sqrt(3)= 2.194 is scale factor
2.2*sqrt(3) # 3.810512
# but Trp-case specificaly yields 3.87 as the average distance between AA (removing first AA outlier)
(sc=3.87/sqrt(3)) # 2.234346 is between 2.2 and 2.3, so use it here in version 2
sc=2.25 # keep things better rounded and not changing with each new protein

pdb$atom=pdb$atom%>%mutate(x=sc*x,y=sc*y,z=sc*z)
pdb$xyz=pdb$xyz*sc
write.pdb(pdb,paste0("trpDia/pdbs/dia",N,".pdb"))


