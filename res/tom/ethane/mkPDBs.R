rm(list=ls())
library(tidyverse) 
library(Rpdb)
(t<-read.pdb("pdbs/ethane.pdb"))
(p<-read.pdb("pdbs/LatDIA.pdb"))
k=c(57,58, 12,22,40, 59,61,63)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/dia.pdb")


(p<-read.pdb("pdbs/LatCUB.pdb"))
k=c(11,14, 2,10,12,13, 23, 15)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/cub.pdb")

(p<-read.pdb("pdbs/LatBCC.pdb"))
k=c(2,15,3,5,9,8,12,14)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/bcc.pdb")

(p<-read.pdb("pdbs/LatFCC.pdb"))
k=c(7,29,6,21,2,28,32,15)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/fcc.pdb")
