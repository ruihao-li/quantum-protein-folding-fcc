rm(list=ls())
library(tidyverse) # for things like mutate()
head(d<-read_delim("trpDia/outs/theta.dat",col_names =c("Theta","RMSD")))
tc=function(sz) theme_classic(base_size=sz);
d%>%ggplot(aes(x=Theta,y=RMSD))+geom_line()+tc(15)
#ggsave("trpDia/outs/rmsd.png",height=3,width=3)
