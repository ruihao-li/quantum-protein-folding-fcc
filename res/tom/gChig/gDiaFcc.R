system("pwd")
library(Peptides) # for plotXVG() views of gmx output
library(NGLVieweR) # for nglviewer looks at gmx inputs (pdbs)
library(dplyr)
base="~/GH/protein-folding-qc/res/tom/gChig" 
setwd(base)
# (old_path <- Sys.getenv("PATH"))
old_path="/home/radivot/miniconda3/bin:/home/radivot/.local/bin:/home/radivot/soft/amber22/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/usr/lib/rstudio/resources/app/bin/quarto/bin:/usr/lib/rstudio/resources/app/bin/postback"
Sys.setenv(PATH = paste("/usr/local/gromacs/bin",old_path, sep = ":"))
system("echo $PATH")
Sys.setenv(GMX_MAXBACKUP=-1)
Sys.setenv(GMX_NO_QUOTES=1)

ks=c(2,3,4,5,6,7,8,9)
(latS=c(paste0("DIA",ks),paste0("FCC",ks)))
dput(latS)
latS=c("DIA3", "DIA4", "DIA5", "DIA6", "DIA7", "DIA8", "DIA9") #,
# "FCC2", "FCC3", "FCC4", "FCC5", "FCC6", "FCC7", "FCC8", "FCC9")
# lat="DIA2"


for (lat in latS) {
  (baseIC=paste0(base,"/",lat,"_350")) 
  system(paste0("mkdir ",baseIC))
  setwd(baseIC) 
  system("mkdir bins")  
  system("mkdir gros")  
  system("mkdir logs")  
  system("mkdir mdps")  
  system("mkdir pdbs")  
  system("mkdir xvgs")  
  system(paste0("gmx pdb2gmx -ignh -f ../../pulchra/pdbChig/",lat,".pdb -ff oplsaa -o gros/pro.gro -p topol.top -i posre.itp -water spce")) #adds H's to pulchra structure
  system("gmx editconf -f gros/pro.gro   -o gros/nbx.gro -c -d 1.0 -bt cubic") # go with making fewer files by treating these gros as temporary
  system("gmx solvate -cp gros/nbx.gro -cs spc216.gro -o gros/sol.gro -p topol.top")
  system("gmx grompp -f ../inputs/ions.mdp -c gros/sol.gro -o bins/ions.tpr -po mdps/mdout.mdp  -p topol.top")
  system("echo 'SOL\n' | gmx genion -s bins/ions.tpr  -o gros/ion.gro -conc 0.15  -pname NA -nname CL -neutral -p topol.top ")
  system("gmx grompp -f ../inputs/minim.mdp -c gros/ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top")
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
  # system(paste0("gmx grompp -f ../inputs/minim.mdp -c gros/",lat,"_Ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top"))
  # system("gmx grompp -f ../inputs/minim.mdp -c gros/ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top")
  system("gmx grompp -f ../inputs/nvt.mdp -c gros/em.gro -r gros/em.gro -p topol.top -o bins/nvt.tpr -po mdps/mdout.mdp") 
  system("gmx mdrun -deffnm bins/nvt")
  system("mv bins/nvt.log logs/nvt.log") 
  system("mv bins/nvt.gro gros/nvt.gro") 
  system("echo 'Temperature' |gmx energy -f bins/nvt.edr -o xvgs/temperature.xvg") 
  plotXVG("xvgs/temperature.xvg")
  
  ######## Raise Water Pressure to Atmospheric ###########
  system("gmx grompp -f ../inputs/npt.mdp -c gros/nvt.gro -r gros/nvt.gro -t bins/nvt.cpt -p topol.top -o bins/npt.tpr -po mdps/mdout.mdp") 
  system("gmx mdrun -deffnm bins/npt")
  system("mv bins/npt.log logs/npt.log") 
  system("mv bins/npt.gro gros/npt.gro") 
  system("echo 'Pressure' | gmx energy -f bins/npt.edr -o xvgs/pressure.xvg") 
  plotXVG("xvgs/pressure.xvg")
  system("echo 'Density' | gmx energy -f bins/npt.edr -o xvgs/density.xvg") 
  plotXVG("xvgs/density.xvg")
  
  ####### now for the production verion. Remove constraints => get rid of -r option
  ######## Hold Room Temp and Pressure and let the protein change shape ###########
  system("gmx grompp -f ../inputs/md.mdp -c gros/npt.gro -t bins/npt.cpt -p topol.top -o bins/md.tpr -po mdps/mdout.mdp") 
  system("gmx mdrun -deffnm bins/md") #uses all 10 cpus by default (30 min on laptop), use -ntmpi 9 to use only 9
  
  #### save these by name to have for analyses later
  system(paste0("mv bins/md.log logs/",lat,".log")) 
  system(paste0("mv bins/md.gro gros/",lat,".gro")) 
  system(paste0("mv bins/md.edr bins/",lat,".edr")) 
  system(paste0("mv bins/md.cpt bins/",lat,".cpt")) 
  system(paste0("mv bins/md.xtc bins/",lat,".xtc")) 
  system(paste0("mv bins/md.tpr bins/",lat,".tpr")) 
  system(paste0("mv bins/em.tpr bins/",lat,"_em.tpr")) 
  system(paste0("mv gros/em.gro gros/",lat,"_em.gro")) 
  
  system(paste0("gmx report-methods -s bins/",lat,".tpr -m logs/",lat)) #makes latex file *.tex
  system(paste0("gmx report-methods -s bins/",lat,"_em.tpr -m logs/",lat,"_em")) #makes latex file *.tex
  
  
  # now make a final pdb that will be compared to the TRU structure
  system(paste0("gmx editconf -f gros/",lat,".gro -o pdbs/tmp0.pdb"))
  system(paste0("grep -v SOL pdbs/tmp0.pdb > pdbs/tmp1.pdb"))
  system(paste0("grep -v NA pdbs/tmp1.pdb  > pdbs/tmp2.pdb"))
  system(paste0("grep -v CL pdbs/tmp2.pdb  > pdbs/",lat,".pdb"))
  system("rm pdbs/tmp*")
  ###### end main computational for loop over lattices and tripeptide center AA numbers 2 to 19
  
  ################### Analyses ##############
  # Select group for centering = 1 (protein), Select group for output = 0 (system)
  # system("gmx trjconv -s bins/md.tpr -f bins/md.xtc -o bins/md_noPBC.xtc -pbc mol -center < inputs/1n0") 
  # system("printf '1\n0\n' | gmx trjconv -s bins/DIA2.tpr -f bins/DIA2.xtc -o bins/DIA2cen.xtc -pbc mol -center") 
  # !printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol # from new py version of tutorial
  system(paste0("printf '1\n0\n' | gmx trjconv -s bins/",lat,".tpr -f bins/",lat,".xtc -o bins/",lat,"_cen.xtc -pbc mol -center")) 
  
  # Select group for least squares fit = 4 (backbone)
  system(paste0("printf '4\n4\n' | gmx rms -s bins/",lat,".tpr -f bins/",lat,"_cen.xtc -o xvgs/rmsd.xvg -tu ns")) 
  plotXVG("xvgs/rmsd.xvg")
  system(paste0("printf '4\n4\n' | gmx rms -s bins/",lat,"_em.tpr -f bins/",lat,"_cen.xtc -o xvgs/rmsd_xtal.xvg -tu ns")) 
  plotXVG("xvgs/rmsd_xtal.xvg")
  system(paste0("printf '1\n' | gmx gyrate -s bins/",lat,".tpr -f bins/",lat,"_cen.xtc -o xvgs/gyrate.xvg")) 
  plotXVG("xvgs/gyrate.xvg")
  
  system(paste0("printf '1\n' | gmx mindist -s bins/",lat,".tpr -f bins/",lat,"_cen.xtc -pi -od xvgs/mindist.xvg")) #makes latex file *.tex
  plotXVG("xvgs/mindist.xvg")
  
  system(paste0('gmx distance  -f bins/',lat,'.xtc -select "atomnr 1,166" -oall xvgs/dists.xvg')) 
  plotXVG("xvgs/dists.xvg")
}

setwd(base)
# NGLVieweR("EXT/pdbs/EXT.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("TRU/pdbs/TRU.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("BCC2/pdbs/BCC2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("BCC10/pdbs/BCC10.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("DIA2/pdbs/DIA2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("FCC5/pdbs/FCC5.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("BCC2/pdbs/BCC250.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("BCC10/pdbs/BCC1050.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("DIA2/pdbs/DIA250.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("FCC5/pdbs/FCC550.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")


