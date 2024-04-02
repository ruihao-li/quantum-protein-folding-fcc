rm(list=ls())
library(tidyverse) 
library(Rpdb)
(t<-read.pdb("pdbs/5AWL_CA.pdb"))
t$atoms$recname="HETATM"
t$atoms$eleid=1:10
t$atoms$elename=" CA"  # need space in front else left-justified => pulchra segmentation fault 
# t$atoms$temp=0.00
t$atoms$occ=1.00
# t$atoms$atomName="C" # try to add element name at end doesn't help
t$conect <- conect(t)

(p<-read.pdb("../TC5b/pdbs/DIA/LatDIA.pdb"))
(d<-read_csv("pdbs/DIA/DIAfits.txt",col_names=F))
for (i in 2:9) {
   # i=2
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("pdbs/DIA/dia",i,".pdb"))
  # write.pdb(t,paste0("pdbs/DIA/fits/dia",i,"test.pdb"))
}

(p<-read.pdb("../TC5b/pdbs/CUB/LatCUB.pdb"))
(d<-read_csv("pdbs/CUB/CUBfits.txt",col_names=F))
for (i in 2:9) {
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("pdbs/CUB/cub",i,".pdb"))
}

(p<-read.pdb("../TC5b/pdbs/BCC/LatBCC.pdb"))
(d<-read_csv("pdbs/BCC/BCCfits.txt",col_names=F))
for (i in 2:9) {
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("pdbs/BCC/bcc",i,".pdb"))
}

(p<-read.pdb("../TC5b/pdbs/FCC/LatFCC.pdb"))
(d<-read_csv("pdbs/FCC/FCCfits.txt",col_names=F))
for (i in 2:9) {
  xyz=p$atoms[as.numeric(d[i-1,]),]%>%select(x1:x3)
  t$atoms[c("x1","x2","x3")]<-xyz
  write.pdb(t,paste0("pdbs/FCC/fcc",i,".pdb"))
}
