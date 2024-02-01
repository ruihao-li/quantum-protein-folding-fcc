rm(list=ls())
library(tidyverse) # for things like mutate()
head(d<-read_delim("trpDia/outs/rmsd6_AB.dat",col_names =c("Index","theta","RMSD")))
tc=function(sz) theme_classic(base_size=sz);
(D=d%>%group_by(Index)%>%summarize(RMSD=min(RMSD)))
D%>%ggplot(aes(x=Index,y=RMSD))+geom_line()+tc(15)+xlab("Index of first AA of pair")
ggsave("trpDia/outs/rmsd6_AB.png",height=3,width=3)

sbb=theme(strip.background=element_blank())
d%>%ggplot(aes(x=theta,y=RMSD))+geom_line()+tc(14)+#ggtitle("Subplot indices are those of first AA of pair aligned")+
  facet_wrap(~Index)+sbb+xlab("Rotation Angle About <1,1,1>")
ggsave("trpDia/outs/rmsd6_facet_AB.png",height=6,width=7.5)

