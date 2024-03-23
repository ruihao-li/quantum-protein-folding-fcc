system("pwd")
# library("here")
# setwd(here())
bccS=c(3,4,8,9,10,11,12)    
lat="BCC"
# i=2
for (i in 3:length(bccS)) {
  (ic=paste0(lat,bccS[i]))
  base="/home/radivot/GH/protein-folding-qc/res/tom/amber" 
  setwd(base)
  # if (0){  # stuff to do just once, to initialize the folder
  #   system("mkdir inputs")
  #   system("cp DIA2/*.in inputs")  #warning, this may wipe out edits made to leap.in ... only do this block on a new folder
  #   system("cp DIA2/*.trajin inputs") 
  #   system("cp DIA2/*.cpptraj inputs") 
  # }
  (baseIC=paste0(base,"/",ic))
  system(paste0("mkdir ",baseIC))
  setwd(baseIC) 
  Sys.setenv(AMBERHOME= "/home/radivot/soft/amber22")
  system("echo $AMBERHOME")
  system("echo $PATH")
  
  con=file("leap.in")
  s1="source leaprc.protein.ff19SB"
  (s2=paste0("TC5b = loadpdb ../../pulchra/pdbTC5b/",ic,".pdb"))
  s3="saveamberparm TC5b TC5b.prmtop TC5b.rst7"
  s4="quit"
  writeLines(c(s1,s2,s3,s4),con)
  system("$AMBERHOME/bin/tleap -s -f leap.in > leap.out")
  system("mpirun -np 12 $AMBERHOME/bin/sander.MPI -O -i ../inputs/min1.in -o min1.out -p TC5b.prmtop -c TC5b.rst7 -r min1.ncrst")
  # Input Files: TC5b.prmtop, TC5b.rst7, min1.in        Output Files: min1.out, min1.ncrst
  system.time(system("mpirun -np 12 $AMBERHOME/bin/sander.MPI -O -i ../inputs/heat1.in -p TC5b.prmtop -c min1.ncrst -r heat1.ncrst -o heat1.out -x heat1.nc"))
  for (i in 2:7)
    print(system.time(system(paste0("mpirun -np 12 $AMBERHOME/bin/pmemd.MPI -O -i ../inputs/heat",i,".in -p TC5b.prmtop -c heat",i-1,".ncrst -r heat",i,".ncrst -o heat",i,".out -x heat",i,".nc"))))
  
  system.time(system("mpirun -np 12 $AMBERHOME/bin/pmemd.MPI -O -i ../inputs/equil.in -p TC5b.prmtop -c heat7.ncrst -r equil1.ncrst -o equil1.out -x equil1.nc")) 
  for (i in 2:10)
    print(system.time(system(paste0("mpirun -np 12 $AMBERHOME/bin/pmemd.MPI -O -i ../inputs/equil.in -p TC5b.prmtop -c equil",i-1,".ncrst -r equil",i,".ncrst -o equil",i,".out -x equil",i,".nc"))))
  
  system("mkdir analysis")
  setwd("analysis")
  system("../../process_mdout.perl ../heat1.out ../heat2.out ../heat3.out ../heat4.out ../heat5.out    ../heat6.out ../heat7.out ../equil1.out ../equil2.out ../equil3.out ../equil4.out ../equil5.out ../equil6.out ../equil7.out ../equil8.out ../equil9.out ../equil10.out")  
  library(readr)
  head(d<-read_table("summary.TEMP",col_names = c("t","temp"))) 
  library(ggplot2)
  library(tidyr)
  # d%>%ggplot(aes(x=t,y=temp))+geom_line()
  # d%>%ggplot(aes(x=t,y=temp))+geom_line()+coord_cartesian(xlim=c(0,50))
  
  head(d<-read_table("summary.EPTOT",col_names = c("t","E"))) 
  # d%>%ggplot(aes(x=t,y=E))+geom_line()
  library(dplyr)
  d=d%>%filter(t>50)
  d%>%ggplot(aes(x=t,y=E))+geom_line()
  (tm=d%>%filter(E==min(E))) 
  system(paste0("grep ' ",tm$t[[1]],"' ../equil*.out > Emin.out"))
  
  con=file("Emin.out")
  (s=readLines(con))
  (nums=as.numeric(stringr::str_extract_all(s, "\\d+")[[1]]))
  (fileNum=nums[1])
  (frameNum=nums[2]/500)
  setwd(baseIC) 
  
  con=file("extract_frame.trajin")
  s1=paste0("trajin equil",fileNum,".nc ",frameNum," ",frameNum)
  s2="trajout lowest_energy_struct.pdb pdb"
  writeLines(c(s1,s2),con)
  
  system("$AMBERHOME/bin/cpptraj TC5b.prmtop < extract_frame.trajin > extract_frame.out")
  
  system("vmd lowest_energy_struct.pdb")  
  
  system("$AMBERHOME/bin/cpptraj TC5b.prmtop < ../inputs/rmsd.trajin > rmsd.out")
  
  head(d<-read_table("rmsd_to_lowest_energy_struct.dat")) 
  names(d) <-c("t","rmsd")
  
  d%>%ggplot(aes(x=t,y=rmsd))+geom_line()
  ggsave("rmsd.pdf")
}
