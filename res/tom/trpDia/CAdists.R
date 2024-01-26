rm(list=ls())
library(tidyverse) # for things like mutate()
library(bio3d) # to read and write pdbs
pdbT <- read.pdb("trpDia/pdbs/2jof_only_alphacarbon.pdb") 
(di=pdbT$atom%>%select(x:z))
L=NULL
for (i in 1:19) L[i]=sqrt((di[i,1]-di[i+1,1])^2+(di[i,2]-di[i+1,2])^2+(di[i,3]-di[i+1,3])^2)
L
mean(L) # 3.864744
sd(L)
# plot(L)
tc=function(sz) theme_classic(base_size=sz);
tibble(i=1:19,L=L)%>%ggplot(aes(x=i,y=L))+geom_line()+xlab("Index of first AA of pair")+
  ylab(expression(paste("Distance between C",alpha," (",ring(A),")" )))+tc(15)
ggsave("trpDia/outs/alphaDists.png",height=3,width=3)

mean(L[-1]) # 3.868206 => 3.87 (first one looks like an outlier so take it out)
mean(L[-(1:9)]) # 3.871271 (remove the whole alpha helix)

#Guo ZY, Kraka E, Cremer D. Description of Local and Global Shape Properties of Protein Helices. J Mol Model. 2013;19(7):2901–2911.
sc=2.3 # The radius, pitch and PPT of an α-helix are 2.3, 5.5 and 3.6 
2.3*sqrt(3) # 3.983717

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3892923/
# distance between consecutive C α atoms are distributed normally with a mean of 3.8 Å and standard deviation of 0.04 Å.
sc=2.194 # 3.8 A in AA space is sqrt(3)=1.7 in diamond space, so 3.8/sqrt(3)= 2.194 is scale facto
2.2*sqrt(3) # 3.810512

3.87/sqrt(3) # 2.234346 is between 2.2 and 2.3, so use it for now

2*sqrt(2)


A=NULL
for (i in 1:18) {
  # i=1
  (v1=t(di[i,]-di[i+1,]))
  (v2=t(t(di[i+2,]-di[i+1,])))
  (v1v2=v2%*%v1)
  sqrt(t(v1)%*%v1)
  (nv1=norm(v1,type="F"))
  (nv2=norm(v2,type="F"))
  A[i]=(180/pi)*acos(v1v2/(nv1*nv2))
}
A
tibble(i=2:19,A=A)%>%ggplot(aes(x=i,y=A))+geom_line()+xlab("Index of 2nd AA of triplet")+
  ylab(expression(paste("C",alpha,"-C",alpha,"-C",alpha," bond angles" )))+tc(15)+geom_hline(yintercept = 109.5)
ggsave("trpDia/outs/angles.png",height=3,width=3)
