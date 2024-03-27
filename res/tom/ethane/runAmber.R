# The barrier to rotation in ethane is about 3.0 kcal/mol and the staggered
# conformation is the lowest-energy conformation. The term used to refer to this
# barrier to rotation is torsional strain.

library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
system("pwd")
Sys.setenv(AMBERHOME= "/home/radivot/soft/amber22")
system("echo $AMBERHOME")
system("echo $PATH")
base="/home/radivot/GH/protein-folding-qc/res/tom/ethane" 
setwd(base)
# system("$AMBERHOME/bin/xleap") #file does not exist, so can't exactly follow https://pharmacy.hsc.wvu.edu/media/1333/ccmmlab3.pdf 
# try doing it via tleap. skipping ahead to antechamber stuff on page 7
# (baseIC=paste0(base,"/AC")) #make an AnteChamber folder
# system(paste0("mkdir ",baseIC))
# setwd(baseIC) 
# system("$AMBERHOME/bin/antechamber -i ../pdbs/ethane.pdb -fi pdb -o et.prepi -fo prepi -c bcc -s 2 -at amber -rn ANE") #makes et.prepi with capital letters

latS=c("TRU","DIA","CUB","BCC","FCC")
lat="CUB"
lat="DIA"
# lat="TRU"
for (lat in latS) {
  (baseIC=paste0(base,"/",lat)) #make an AnteCham folder using gaff2 (lower case letters in prepi file)
  system(paste0("mkdir ",baseIC))
  setwd(baseIC) 
  # system(paste0("$AMBERHOME/bin/pdb4amber -i ../pdbs/",tolower(lat),".pdb -o ",tolower(lat),"2.pdb"))
  # system(paste0("$AMBERHOME/bin/antechamber -j 5 -at sybyl -dr no -i ../pdbs/",tolower(lat),".pdb -fi pdb -o ../pdbs/",tolower(lat),"2.pdb -fo pdb"))
  system(paste0("$AMBERHOME/bin/antechamber -i ../pdbs/",tolower(lat),".pdb -fi pdb -o et.prepi -fo prepi -c bcc -s 2 -at gaff2 -rn ANE"))
  system("rm A*")
  system("rm P*")
  system("rm sqm*")  #check sqm.out first, looks fine, so delete
  # system("$AMBERHOME/bin/parmchk2 -i et.prepi -f prepi -o et.frcmod") #regular parmchk not found. All empty good, no probs
  
  con=file("leap.in")
  s1="source leaprc.gaff2"
  (s2="loadamberprep et.prepi")
  # (s2=paste0("Et = loadpdb ../pdbs/",tolower(lat),".pdb"))
  # (s3="loadamberparams et.frcmod") # not needed since empty (all were good)
  # s3="ANE = loadpdb ../pdbs/ethane.pdb"  # still fails with UNK
  s3="et = loadpdb NEWPDB.PDB" 
  # s4="et = copy ANE"
  # s5="saveoff et et.lib"
  # s6="savepdb et et.pdb"
  # saveamberparm tbutane tbutane.top tbutane.crd
  s7="saveamberparm et et.prmtop et.rst7"
  s8="quit"
  writeLines(c(s1,s2,s3,s7,s8),con)
  # writeLines(c(s1,s2,s3,s4,s5,s6,s7,s8),con)
  close(con)
  system("$AMBERHOME/bin/tleap -s -f leap.in > leap.out")
  system("$AMBERHOME/bin/sander.MPI -O -i ../../amber/inputs/min1.in -o min1.out -p et.prmtop -c et.rst7 -r min1.ncrst")
  system.time(system("$AMBERHOME/bin/sander.MPI -O -i ../../amber/inputs/heat1.in -p et.prmtop -c min1.ncrst -r heat1.ncrst -o heat1.out -x heat1.nc"))
  for (i in 2:7)
    print(system.time(system(paste0("mpirun -np 1 $AMBERHOME/bin/sander.MPI -O -i ../../amber/inputs/heat",i,".in -p et.prmtop -c heat",i-1,".ncrst -r heat",i,".ncrst -o heat",i,".out -x heat",i,".nc"))))
  
  system.time(system("mpirun -np 1 $AMBERHOME/bin/sander.MPI -O -i ../../amber/inputs/equil.in -p et.prmtop -c heat7.ncrst -r equil1.ncrst -o equil1.out -x equil1.nc")) #10 secs
  for (i in 2:10)
    print(system.time(system(paste0("mpirun -np 1 $AMBERHOME/bin/sander.MPI -O -i ../../amber/inputs/equil.in -p et.prmtop -c equil",i-1,".ncrst -r equil",i,".ncrst -o equil",i,".out -x equil",i,".nc"))))
  
  system("rm N*")
  
  system("mkdir analysis")
  setwd("analysis")
  system("../../../amber/process_mdout.perl ../heat1.out ../heat2.out ../heat3.out ../heat4.out ../heat5.out    ../heat6.out ../heat7.out ../equil1.out ../equil2.out ../equil3.out ../equil4.out ../equil5.out ../equil6.out ../equil7.out ../equil8.out ../equil9.out ../equil10.out")
  # system("../../../amber/process_mdout.perl ../heat1.out ../heat2.out ../heat3.out ../heat4.out ../heat5.out ../heat6.out ../heat7.out ../equil1.out ../equil2.out ../equil3.out")  
  head(d<-read_table("summary.EPTOT",col_names = c("t","E"))) 
  d=d%>%filter(t>50)
  (tm=d%>%filter(E==min(E))) 
  system(paste0("grep ' ",tm$t[[1]],"' ../equil*.out > Emin.out"))
  
  con=file("Emin.out")
  (s=readLines(con))
  close(con)
  (nums=as.numeric(stringr::str_extract_all(s, "\\d+")[[1]]))
  (fileNum=nums[1])
  (frameNum=nums[2]/500)
  setwd(baseIC) 
  
  con=file("extract_frame.trajin")
  s1=paste0("trajin equil",fileNum,".nc ",frameNum," ",frameNum)
  s2="trajout lowest_energy_struct.pdb pdb"
  writeLines(c(s1,s2),con)
  close(con)
  system("$AMBERHOME/bin/cpptraj et.prmtop < extract_frame.trajin > extract_frame.out")
  system("$AMBERHOME/bin/cpptraj et.prmtop < ../../amber/inputs/rmsdCH.trajin > rmsd.out")
  
  head(d<-read_table("rmsd_to_tru.dat")) 
  names(d) <-c("t","rmsd")
  
  d%>%filter(t<500)%>%ggplot(aes(x=t,y=rmsd))+geom_line()
  ggsave("rmsd.pdf")
  # system("rm lowest_energy_struct.pdb")
  system("rm rmsd_to_tru.dat")
}




