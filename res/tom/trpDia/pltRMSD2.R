rm(list=ls())
library(tidyverse) # for things like mutate()
head(d<-read_delim("task1/outs/rmsd2.dat",col_names =c("Index","RMSD")))
tc=function(sz) theme_classic(base_size=sz);
d%>%ggplot(aes(x=Index,y=RMSD))+geom_line()+tc(15)+xlab("Index of first AA of pair")
ggsave("task1/outs/rmsd2.png",height=3,width=3)
