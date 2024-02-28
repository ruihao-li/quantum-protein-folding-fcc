# get pulchra from  https://github.com/euplotes/pulchra
# and compile it like this
# system("cc -O3 -o pulchra pulchra.c pulchra_data.c -lm")
# to get pulchra to run from anywhere, soft link it into a folder in the path, e.g. 
# cd /usr/local/bin
# sudo ln -s /Users/radivot/ccf/quantum/pulchra/pulchra pulchra

for (i in c("DIA","CUB","BCC","FCC")) 
  for (j in 2:19) {
    system(paste0("pulchra -v trpLat/pdbs/",i,"/fits/",tolower(i),j,".pdb")) 
    system(paste0("mv trpLat/pdbs/",i,"/fits/",tolower(i),j,".rebuilt.pdb pulchra/pdbs/",i,j,".pdb")) 
  } 
