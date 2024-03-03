library(Peptides) # for plotXVG() views of gmx output
library(NGLVieweR) # for nglviewer looks at gmx inputs (pdbs)
library(tidyverse)
setwd("gmx") # switch to gmx to keep topol.top and posre.top in the gmx folder
Sys.setenv(GMX_MAXBACKUP=-1)
Sys.setenv(GMX_NO_QUOTES=1)
system("set")
NGLVieweR("pdbs/2jof_F1.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("../pulchra/pdbs/DIA2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("../pulchra/pdbs/CUB2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("../pulchra/pdbs/BCC2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("../pulchra/pdbs/FCC2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")

# i=2
# lat="DIA"  #set things up for for-loops later
for (i in 3:19)
  for (lat in c("DIA")) {
#   for (lat in c("FCC")) {
    # system(paste0("gmx pdb2gmx -ignh -f ../pulchra/pdbs/",lat,i,".pdb -ff oplsaa -o gros/",lat,i,"_Pro.gro -p topol.top -i posre.itp -water spce"))
    system(paste0("gmx pdb2gmx -ignh -f ../pulchra/pdbs/",lat,i,".pdb -ff oplsaa -o gros/pro.gro -p topol.top -i posre.itp -water spce"))
    # different options used in https://tutorials.gromacs.org/docs/md-intro-tutorial.html are -ff charm27 and -water tip3p
    # system(paste0("gmx editconf -f gros/",lat,i,"_Pro.gro   -o gros/",lat,i,"_Nbx.gro -c -d 1.0 -bt cubic"))
    system("gmx editconf -f gros/pro.gro   -o gros/nbx.gro -c -d 1.0 -bt cubic") # go with making fewer files by treating these gros as temporary
    # system(paste0("gmx solvate -cp gros/",lat,i,"_Nbx.gro -cs spc216.gro -o gros/",lat,i,"_Sol.gro -p topol.top"))
    system("gmx solvate -cp gros/nbx.gro -cs spc216.gro -o gros/sol.gro -p topol.top")
    # system("touch mdps/ions.mdp") #do it once to make an empty dummy file input for grompp
    # system(paste0("gmx grompp -f mdps/ions.mdp -c gros/",lat,i,"_Sol.gro -o bins/ions.tpr -po mdps/mdout.mdp  -p topol.top"))
    system("gmx grompp -f mdps/ions.mdp -c gros/sol.gro -o bins/ions.tpr -po mdps/mdout.mdp  -p topol.top")
    # system(paste0("echo 'SOL\n' | gmx genion -s bins/ions.tpr  -o gros/",lat,i,"_Ion.gro -conc 0.15  -pname NA -nname CL -neutral -p topol.top "))
    system("echo 'SOL\n' | gmx genion -s bins/ions.tpr  -o gros/ion.gro -conc 0.15  -pname NA -nname CL -neutral -p topol.top ")
    # system(paste0("gmx grompp -f mdps/minim.mdp -c gros/",lat,i,"_Ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top"))
    system("gmx grompp -f mdps/minim.mdp -c gros/ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top")
    system("tail -8 topol.top" )
    
    #### many files in the next 3 chunks will be overwritten in for loops (treated as temporary files)
    ######## Energy Minimization ###########
    system("gmx mdrun -deffnm bins/em") #deffnm = default file name, -v makes it verbose (skip in real runs)
    ## makes all of these files
    # em.log: ASCII-text log file of the EM process
    # em.edr: Binary energy file
    # em.trr: Binary full-precision trajectory
    # em.gro: Energy-minimized structure
    system("mv bins/em.log logs/em.log") 
    system("mv bins/em.gro gros/em.gro") 
    # system("gmx energy -f bins/em.edr -o xvgs/potential.xvg < inputs/10_0")
    system("echo 'Potential\n0\n' | gmx energy -f bins/em.edr -o xvgs/potential.xvg")
    plotXVG("xvgs/potential.xvg") # from R package Peptides
    
    ######## Raise Water Temperature to Room Temp ###########
    system(paste0("gmx grompp -f mdps/minim.mdp -c gros/",lat,i,"_Ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top"))
    system("gmx grompp -f mdps/nvt.mdp -c gros/em.gro -r gros/em.gro -p topol.top -o bins/nvt.tpr -po mdps/mdout.mdp") 
    system("gmx mdrun -deffnm bins/nvt")
    system("mv bins/nvt.log logs/nvt.log") 
    system("mv bins/nvt.gro gros/nvt.gro") 
    system("echo 'Temperature' |gmx energy -f bins/nvt.edr -o xvgs/temperature.xvg") 
    plotXVG("xvgs/temperature.xvg")
    
    ######## Raise Water Pressure to Atmospheric ###########
    system("gmx grompp -f mdps/npt.mdp -c gros/nvt.gro -r gros/nvt.gro -t bins/nvt.cpt -p topol.top -o bins/npt.tpr -po mdps/mdout.mdp") 
    system("gmx mdrun -deffnm bins/npt")
    system("mv bins/npt.log logs/npt.log") 
    system("mv bins/npt.gro gros/npt.gro") 
    system("echo 'Pressure' | gmx energy -f bins/npt.edr -o xvgs/pressure.xvg") 
    plotXVG("xvgs/pressure.xvg")
    system("echo 'Density' | gmx energy -f bins/npt.edr -o xvgs/density.xvg") 
    plotXVG("xvgs/density.xvg")
    
    ####### now for the production verion. Remove constraints => get rid of -r option
    ######## Hold Room Temp and Pressure and let the protein change shape ###########
    system("gmx grompp -f mdps/md.mdp -c gros/npt.gro -t bins/npt.cpt -p topol.top -o bins/md.tpr -po mdps/mdout.mdp") 
    system("gmx mdrun -deffnm bins/md") #uses all 10 cpus by default (30 min on laptop), use -ntmpi 9 to use only 9
    
    #### save these by name to have for analyses later
    system(paste0("mv bins/md.log logs/",lat,i,".log")) 
    system(paste0("mv bins/md.gro gros/",lat,i,".gro")) 
    system(paste0("mv bins/md.edr bins/",lat,i,".edr")) 
    system(paste0("mv bins/md.cpt bins/",lat,i,".cpt")) 
    system(paste0("mv bins/md.xtc bins/",lat,i,".xtc")) 
    system(paste0("mv bins/md.tpr bins/",lat,i,".tpr")) 
    system(paste0("mv bins/em.tpr bins/",lat,i,"_em.tpr")) 
    system(paste0("mv gros/em.gro gros/",lat,i,"_em.gro")) 
    
    system(paste0("gmx report-methods -s bins/",lat,i,".tpr -m logs/",lat,i)) #makes latex file *.tex
    system(paste0("gmx report-methods -s bins/",lat,i,"_em.tpr -m logs/",lat,i,"_em")) #makes latex file *.tex
    
    
    # now make a final pdb that will be compared to the trp-cage xtal structure 2jof_F1.pdb
    system(paste0("gmx editconf -f gros/",lat,i,".gro -o pdbs/tmp0.pdb"))
    system(paste0("grep -v SOL pdbs/tmp0.pdb > pdbs/tmp1.pdb"))
    system(paste0("grep -v NA pdbs/tmp1.pdb  > pdbs/tmp2.pdb"))
    system(paste0("grep -v CL pdbs/tmp2.pdb  > pdbs/",lat,i,".pdb"))
    system("rm pdbs/tmp*")
  }



###### end main computational for loop over lattices and tripeptide center AA numbers 2 to 19

################### Analyses ##############
# Select group for centering = 1 (protein), Select group for output = 0 (system)
# system("gmx trjconv -s bins/md.tpr -f bins/md.xtc -o bins/md_noPBC.xtc -pbc mol -center < inputs/1n0") 
# system("printf '1\n0\n' | gmx trjconv -s bins/DIA2.tpr -f bins/DIA2.xtc -o bins/DIA2cen.xtc -pbc mol -center") 
# !printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol # from new py version of tutorial
system(paste0("printf '1\n0\n' | gmx trjconv -s bins/",lat,i,".tpr -f bins/",lat,i,".xtc -o bins/",lat,i,"_cen.xtc -pbc mol -center")) 

# Select group for least squares fit = 4 (backbone)
system(paste0("printf '4\n4\n' | gmx rms -s bins/",lat,i,".tpr -f bins/",lat,i,"_cen.xtc -o xvgs/rmsd.xvg -tu ns")) 
plotXVG("xvgs/rmsd.xvg")
system(paste0("printf '4\n4\n' | gmx rms -s bins/",lat,i,"_em.tpr -f bins/",lat,i,"_cen.xtc -o xvgs/rmsd_xtal.xvg -tu ns")) 
plotXVG("xvgs/rmsd_xtal.xvg")
system(paste0("printf '1\n' | gmx gyrate -s bins/",lat,i,".tpr -f bins/",lat,i,"_cen.xtc -o xvgs/gyrate.xvg")) 
plotXVG("xvgs/gyrate.xvg")

system(paste0("printf '1\n' | gmx mindist -s bins/",lat,i,".tpr -f bins/",lat,i,"_cen.xtc -pi -od xvgs/mindist.xvg")) #makes latex file *.tex
plotXVG("xvgs/mindist.xvg")

######## Notes (to get Trp 6 distances to prolines 12 and 18) ################
# gmx distance calculates distances between pairs of positions as a function of time. NOT what we want here
#
#########  Option 1 ##########
# gmx mindist computes the distance between one group and a number of other
# groups. Both the minimum distance (between any pair of atoms from the
# respective groups) and the number of contacts within a given distance are
# written to two separate output files. With the -group option a contact of an
# atom in another group with multiple atoms in the first group is counted as one
# contact instead of as multiple contacts. With -or, minimum distances to each
# residue in the first group are determined and plotted as a function of residue
# number.
#
# With option -pi the minimum distance of a group to its periodic image is
# plotted. This is useful for checking if a protein has seen its periodic image
# during a simulation. Only one shift in each direction is considered, giving a
# total of 26 shifts. Note that periodicity information is required from the
# file supplied with with -s, either as a .tpr file or a .pdb file with CRYST1
# fields. It also plots the maximum distance within the group and the lengths of
# the three box vectors.
#
#######  Option 2 ##########

# gmx pairdist calculates pairwise distances between one reference selection
# (given with -ref) and one or more other selections (given with -sel). It can
# calculate either the minimum distance (the default), or the maximum distance
# (with -type max). Distances to each selection provided with -sel are computed
# independently.
#
# By default, the global minimum/maximum distance is computed. To compute more
# distances (e.g., minimum distances to each residue in -ref), use -refgrouping
# and/or -selgrouping to specify how the positions within each selection should
# be grouped.
#
# Computed distances are written to the file specified with -o. If there are N
# groups in -ref and M groups in the first selection in -sel, then the output
# contains N*M columns for the first selection. The columns contain distances
# like this: r1-s1, r2-s1, …, r1-s2, r2-s2, …, where rn is the n’th group in
# -ref and sn is the n’th group in the other selection. The distances for the
# second selection comes as separate columns after the first selection, and so
# on. If some selections are dynamic, only the selected positions are used in
# the computation but the same number of columns is always written out. If there
# are no positions contributing to some group pair, then the cutoff value is
# written (see below).
#
# cutoff sets a cutoff for the computed distances. If the result would contain a
# distance over the cutoff, the cutoff value is written to the output file
# instead. By default, no cutoff is used, but if you are not interested in
# values beyond a cutoff, or if you know that the minimum distance is smaller
# than a cutoff, you should set this option to allow the tool to use grid-based
# searching and be significantly faster.
#
# If you want to compute distances between fixed pairs, gmx distance may be a
# more suitable tool.


# ############# Junk yard below #######
# #### Unfinished start on doing it in python
# import os
# i=2
# lat="DIA"
# os.system("gmx pdb2gmx -ignh -f ../pulchra/pdbs/"+lat+str(i)+".pdb -ff oplsaa -o gros/"+lat+str(i)+"_Pro.gro -p topol.top -i posre.itp -water spce")
# os.system("gmx editconf -f gros/"+lat+str(i)+"_Pro.gro  -o gros/"+lat+str(i)+"_Nbx.gro -c -d 1.0 -bt cubic")
# os.system("gmx solvate -cp gros/"+lat+str(i)+"_Nbx.gro -cs spc216.gro -o gros/"+lat+str(i)+"_Sol.gro -p topol.top")
# # os.system("touch mdps/ions.mdp") #do it once to make an empty dummy file input for grompp
# os.system("gmx grompp -f mdps/ions.mdp -c gros/"+lat+str(i)+"_Sol.gro -o bins/ions.tpr -po mdps/mdout.mdp  -p topol.top")
# os.system("echo \"SOL\n\"  > gmx genion -s bins/ions.tpr  -o gros/"+lat+str(i)+"_Ion.gro -conc 0.15 -pname NA -nname CL -neutral -p topol.top ")
# !tail -8 topol.top #74 SOL replaced  by NA and CL to make the ion concentration right (and total charge neutral)


