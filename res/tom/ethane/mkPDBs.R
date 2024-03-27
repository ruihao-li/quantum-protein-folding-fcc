rm(list=ls())  #Edit output of this by hand to match tru.pdb format
# library(tidyverse)
library(Rpdb)
setwd("~/GH/protein-folding-qc/res/tom/ethane")

(t<-read.pdb("pdbs/tru.pdb"))
t$atoms$elename=paste0(" ",t$atoms$elename)
# " CA"  # need space in front else left-justified => pulchra segmentation fault
t$conect=NULL
s1="CONECT    1    3    4    2    5"
s2="CONECT    3    1"
s3="CONECT    4    1"
s4="CONECT    5    1"
s5="CONECT    2    1    6    7    8"
s6="CONECT    6    2"
s7="CONECT    7    2"
s8="CONECT    8    2"
s9="END   "
con=file("pdbs/end.txt")
writeLines(c(s1,s2,s3,s4,s5,s6,s7,s8,s9),con)
close(con)


(p<-read.pdb("pdbs/LatDIA.pdb"))
k=c(57,58, 12,22,40, 59,61,63)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/dia2.pdb")
system("cat pdbs/dia2.pdb pdbs/end.txt > pdbs/dia.pdb")

(p<-read.pdb("pdbs/LatCUB.pdb"))
k=c(11,14, 2,10,12,13, 23, 15)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/cub2.pdb")
system("cat pdbs/cub2.pdb pdbs/end.txt > pdbs/cub.pdb")

(p<-read.pdb("pdbs/LatBCC.pdb"))
k=c(2,15,3,5,9,8,12,14)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/bcc2.pdb")
system("cat pdbs/bcc2.pdb pdbs/end.txt > pdbs/bcc.pdb")

(p<-read.pdb("pdbs/LatFCC.pdb"))
k=c(7,29,6,21,2,28,32,15)
xyz=p$atoms[k,]%>%select(x1:x3)
t$atoms[c("x1","x2","x3")]<-xyz
write.pdb(t,"pdbs/fcc2.pdb")
system("cat pdbs/fcc2.pdb pdbs/end.txt > pdbs/fcc.pdb")

system("rm pdbs/*2.pdb")
system("rm pdbs/end.txt")