library(Peptides) # for plotXVG() views of gmx output
library(NGLVieweR) # for nglviewer looks at gmx inputs (pdbs)
library(dplyr)
base="~/GH/protein-folding-qc/res/tom/gAlaDiP" 
setwd(base)
# (old_path <- Sys.getenv("PATH"))
old_path="/home/radivot/miniconda3/bin:/home/radivot/.local/bin:/home/radivot/soft/amber22/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/usr/lib/rstudio/resources/app/bin/quarto/bin:/usr/lib/rstudio/resources/app/bin/postback"
Sys.setenv(PATH = paste("/usr/local/gromacs/bin",old_path, sep = ":"))
system("echo $PATH")
Sys.setenv(GMX_MAXBACKUP=-1)
Sys.setenv(GMX_NO_QUOTES=1)
NGLVieweR("inputs/alaDiP.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")

(baseIC=paste0(base,"/350")) 
system(paste0("mkdir ",baseIC))
setwd(baseIC) 
system("mkdir bins")  
system("mkdir gros")  
system("mkdir logs")  
system("mkdir mdps")  
system("mkdir pdbs")  
system("mkdir xvgs")  
######## Energy Minimization ###########
system("gmx mdrun -s ../inputs/topol_tpr -nsteps 5000000 -ntomp 1   -deffnm bins/em") #deffnm = default file name, -v makes it verbose (skip in real runs)
system("mv bins/em.log logs/em.log") 
system("mv bins/em.gro gros/em.gro") 
system("echo 'Potential\n0\n' | gmx energy -f bins/em.edr -o xvgs/potential.xvg")
plotXVG("xvgs/potential.xvg") # from R package Peptides
system("echo 'Temperature' |gmx energy -f bins/em.edr -o xvgs/temperature.xvg") 
plotXVG("xvgs/temperature.xvg")
system("echo 'Pressure' | gmx energy -f bins/em.edr -o xvgs/pressure.xvg") 
plotXVG("xvgs/pressure.xvg")
system("echo 'Density' | gmx energy -f bins/em.edr -o xvgs/density.xvg") 
plotXVG("xvgs/density.xvg")
system("printf '4\n4\n' | gmx rms -s ../inputs/alaDiP.pdb -f bins/em.xtc -o xvgs/rmsd.xvg -tu ns") 
plotXVG("xvgs/rmsd.xvg")
system("printf '1\n' | gmx gyrate -s ../inputs/alaDiP.pdb -f bins/em.xtc -o xvgs/gyrate.xvg") 
plotXVG("xvgs/gyrate.xvg")

