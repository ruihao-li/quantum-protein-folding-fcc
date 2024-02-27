rm(list=ls())
library(tidyverse) 
library(Rpdb)
(t<-read.pdb("trpDia/pdbs/trpCen.pdb"))
t$atoms$recname="HETATM"
t$atoms$eleid=1:20
t$atoms$elename=" CA"  # need space in front else left-justified => pulchra segmentation fault 
t$atoms$temp=0.00
t$atoms$occ=1.00
t$conect <- conect(t)

(p<-read.pdb("trpLat/pdbs/DIA/LatDIA.pdb"))
(d<-read_csv("trpLat/pdbs/DIA/DIAfits.txt",col_names=F))
for (i in 2:19) {
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("trpLat/pdbs/DIA/fits/dia",i,".pdb"))
}

(p<-read.pdb("trpLat/pdbs/CUB/LatCUB.pdb"))
(d<-read_csv("trpLat/pdbs/CUB/CUBfits.txt",col_names=F))
for (i in 2:19) {
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("trpLat/pdbs/CUB/fits/cub",i,".pdb"))
}

(p<-read.pdb("trpLat/pdbs/BCC/LatBCC.pdb"))
(d<-read_csv("trpLat/pdbs/BCC/BCCfits.txt",col_names=F))
for (i in 2:19) {
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("trpLat/pdbs/BCC/fits/bcc",i,".pdb"))
}

(p<-read.pdb("trpLat/pdbs/FCC/LatFCC.pdb"))
(d<-read_csv("trpLat/pdbs/FCC/FCCfits.txt",col_names=F))
for (i in 2:19) {
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("trpLat/pdbs/FCC/fits/fcc",i,".pdb"))
}
