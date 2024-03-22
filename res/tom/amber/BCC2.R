system("pwd")
setwd("/home/radivot/GH/protein-folding-qc/res/tom/amber/BCC2") 
if (0){  # stuff to do just once, to initialize the folder
  system("mkdir analysis")
  system("cp ../DIA2/*.in .")  #warning, this may wipe out edits made to leap.in ... only do this block on a new folder
  system("cp ../DIA2/*.trajin .") 
  system("cp ../DIA2/*.cpptraj .") 
}
Sys.setenv(AMBERHOME= "/home/radivot/soft/amber22")
system("echo $AMBERHOME")
system("echo $PATH")
system("$AMBERHOME/bin/tleap -s -f leapB2.in > leapB2.out")
system("mpirun -np 12 $AMBERHOME/bin/sander.MPI -O -i min1.in -o min1.out -p TC5b.prmtop -c TC5b.rst7 -r min1.ncrst")


# Input Files: TC5b.prmtop, TC5b.rst7, min1.in        Output Files: min1.out, min1.ncrst
# system("$AMBERHOME/bin/ambpdb -p TC5b.prmtop -c min1.ncrst > min1.pdb")
system.time(system("mpirun -np 12 $AMBERHOME/bin/sander.MPI -O -i heat1.in -p TC5b.prmtop -c min1.ncrst -r heat1.ncrst -o heat1.out -x heat1.nc"))
for (i in 2:7)
  print(system.time(system(paste0("mpirun -np 12 $AMBERHOME/bin/pmemd.MPI -O -i heat",i,".in -p TC5b.prmtop -c heat",i-1,".ncrst -r heat",i,".ncrst -o heat",i,".out -x heat",i,".nc"))))

system.time(system("mpirun -np 12 $AMBERHOME/bin/pmemd.MPI -O -i equil.in -p TC5b.prmtop -c heat7.ncrst -r equil1.ncrst -o equil1.out -x equil1.nc"))  #works
for (i in 2:10)
  print(system.time(system(paste0("mpirun -np 12 $AMBERHOME/bin/pmemd.MPI -O -i equil.in -p TC5b.prmtop -c equil",i-1,".ncrst -r equil",i,".ncrst -o equil",i,".out -x equil",i,".nc"))))
# system("vmd -parm7 TC5b.prmtop -netcdf heat10.nc")  

setwd("analysis")
system("../../process_mdout.perl ../heat1.out ../heat2.out ../heat3.out ../heat4.out ../heat5.out    ../heat6.out ../heat7.out ../equil1.out ../equil2.out ../equil3.out ../equil4.out ../equil5.out ../equil6.out ../equil7.out ../equil8.out ../equil9.out ../equil10.out")  
library(readr)
head(d<-read_table("summary.TEMP",col_names = c("t","temp"))) 
library(ggplot2)
library(tidyr)
d%>%ggplot(aes(x=t,y=temp))+geom_line()
d%>%ggplot(aes(x=t,y=temp))+geom_line()+coord_cartesian(xlim=c(0,50))

head(d<-read_table("summary.EPTOT",col_names = c("t","E"))) 
d%>%ggplot(aes(x=t,y=E))+geom_line()
library(dplyr)
d=d%>%filter(t>50)
d%>%ggplot(aes(x=t,y=E))+geom_line()
(tm=d%>%filter(E==min(E))) 
(s=paste0("grep ",tm$t[[1]]," ../*.out"))
system(s)
system("grep 30892 ../*.out")
421000/500 # frame 842 in file 7

setwd("/home/radivot/GH/protein-folding-qc/res/tom/amber/BCC2") 
system("$AMBERHOME/bin/cpptraj TC5b.prmtop < extract_frame7.trajin > extract_frame7.out")

system("vmd lowest_energy_struct.pdb")  #looks good

system("$AMBERHOME/bin/cpptraj TC5b.prmtop < rmsd.trajin > rmsd.out")

head(d<-read_table("rmsd_to_lowest_energy_struct.dat")) 
names(d) <-c("t","rmsd")

d%>%ggplot(aes(x=t,y=rmsd))+geom_line()
# ggsave("rmsd.pdf")

system("$AMBERHOME/bin/cpptraj TC5b.prmtop <nc_to_binpos.cpptraj >nc_to_binpos.out")
system("$AMBERHOME/bin/cpptraj TC5b.prmtop <average_31-32_crd.cpptraj")
system("$AMBERHOME/bin/cpptraj TC5b.prmtop <measure_trp_angles.cpptraj")

system("vmd equil_average_31-32.pdb")  ## same

