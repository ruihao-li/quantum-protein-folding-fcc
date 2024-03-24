library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
setwd("~/GH/protein-folding-qc/res/tom/amber/trp19") 
head(d<-read_table("rmsd_to_lowest_energy_struct.dat")) 
names(d) <-c("t","rmsd")
d%>%ggplot(aes(x=t/1000,y=rmsd))+geom_line()+theme_classic()+xlab("Time (ns)")+ylab("RMSD")
ggsave("rmsd.png",width=6,height=3) #lots of data => pdfs are big

setwd("~/GH/protein-folding-qc/res/tom/amber/BCC2") 
head(d<-read_table("rmsd_to_lowest_energy_struct.dat")) 
names(d) <-c("t","rmsd")
d%>%ggplot(aes(x=t/1000,y=rmsd))+geom_line()+theme_classic()+xlab("Time (ns)")+ylab("RMSD")
ggsave("rmsd.png",width=6,height=3) #lots of data => pdfs are big

setwd("~/GH/protein-folding-qc/res/tom/amber/BCC10") 
head(d<-read_table("rmsd_to_lowest_energy_struct.dat")) 
names(d) <-c("t","rmsd")
d%>%ggplot(aes(x=t/1000,y=rmsd))+geom_line()+theme_classic()+xlab("Time (ns)")+ylab("RMSD")
ggsave("rmsd.png",width=6,height=3) #lots of data => pdfs are big


# this worked in ubuntu but not on Mac, do from OS
# system("vmd lowest_energy_struct.pdb")  ### excellent fit of Trp in Cage!!!
