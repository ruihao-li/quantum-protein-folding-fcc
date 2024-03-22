system("pwd")
setwd("/home/radivot/GH/protein-folding-qc/res/tom/amber/FCC10") 
Sys.setenv(AMBERHOME= "/home/radivot/soft/amber22")
system("echo $AMBERHOME")
system("echo $PATH")
system("$AMBERHOME/bin/tleap -s -f leapF10.in > leapF10.out")
# FCC10 has the biggest RMSD (of FCCs), and it too fails, just like FCC5, which has the smallest RMSD 
system("mpirun -np 12 $AMBERHOME/bin/sander.MPI -O -i min1.in -o min1.out -p TC5b.prmtop -c TC5b.rst7 -r min1.ncrst")
