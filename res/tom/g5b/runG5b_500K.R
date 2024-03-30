library(Peptides) # for plotXVG() views of gmx output
library(NGLVieweR) # for nglviewer looks at gmx inputs (pdbs)
library(dplyr)
base="~/GH/protein-folding-qc/res/tom/g5b" 
setwd(base)
# (old_path <- Sys.getenv("PATH"))
old_path="/home/radivot/miniconda3/bin:/home/radivot/.local/bin:/home/radivot/soft/amber22/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/usr/lib/rstudio/resources/app/bin/quarto/bin:/usr/lib/rstudio/resources/app/bin/postback"
Sys.setenv(PATH = paste("/usr/local/gromacs/bin",old_path, sep = ":"))
system("echo $PATH")
Sys.setenv(GMX_MAXBACKUP=-1)
Sys.setenv(GMX_NO_QUOTES=1)
latS="TRU"  # 50 ns =>2 hours   
latS="EXT"  # 50 ns => 12+ hours due to bigger box of water (in top see 18529 waters vs 2662 in TRU)
latS=c("BCC2","BCC10","DIA2","FCC5") # based on amber, only last two should fail (go back to 1 ns first just to test this, no go back up to 50 ns)

for (lat in latS) {
  (baseIC=paste0(base,"/",lat,"_500")) 
  system(paste0("mkdir ",baseIC))
  setwd(baseIC) 
  system("mkdir bins")  
  system("mkdir gros")  
  system("mkdir logs")  
  system("mkdir mdps")  
  system("mkdir pdbs")  
  system("mkdir xvgs")  
  system(paste0("gmx pdb2gmx -ignh -f ../../pulchra/pdbTC5b/",lat,".pdb -ff oplsaa -o gros/pro.gro -p topol.top -i posre.itp -water spce"))
  system("gmx editconf -f gros/pro.gro   -o gros/nbx.gro -c -d 1.0 -bt cubic") # go with making fewer files by treating these gros as temporary
  system("gmx solvate -cp gros/nbx.gro -cs spc216.gro -o gros/sol.gro -p topol.top")
  system("gmx grompp -f ../inputs/ions.mdp -c gros/sol.gro -o bins/ions.tpr -po mdps/mdout.mdp  -p topol.top")
  system("echo 'SOL\n' | gmx genion -s bins/ions.tpr  -o gros/ion.gro -conc 0.15  -pname NA -nname CL -neutral -p topol.top ")
  system("gmx grompp -f ../inputs/minim.mdp -c gros/ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top")
  system("tail -8 topol.top" )
  ######## Energy Minimization ###########
  system("gmx mdrun -deffnm bins/em") #deffnm = default file name, -v makes it verbose (skip in real runs)
  system("mv bins/em.log logs/em.log") 
  system("mv bins/em.gro gros/em.gro") 
  system("echo 'Potential\n0\n' | gmx energy -f bins/em.edr -o xvgs/potential.xvg")
  plotXVG("xvgs/potential.xvg") # from R package Peptides
  
  ######## Raise Water Temperature to Room Temp ###########
  system(paste0("gmx grompp -f ../inputs/minim.mdp -c gros/",lat,"_Ion.gro -o bins/em.tpr -po mdps/mdout.mdp  -p topol.top"))
  system("gmx grompp -f ../inputs/nvt500.mdp -c gros/em.gro -r gros/em.gro -p topol.top -o bins/nvt.tpr -po mdps/mdout.mdp") 
  system("gmx mdrun -deffnm bins/nvt")
  system("mv bins/nvt.log logs/nvt.log") 
  system("mv bins/nvt.gro gros/nvt.gro") 
  system("echo 'Temperature' |gmx energy -f bins/nvt.edr -o xvgs/temperature.xvg") 
  plotXVG("xvgs/temperature.xvg")
  
  ######## Raise Water Pressure to Atmospheric ###########
  system("gmx grompp -f ../inputs/npt500.mdp -c gros/nvt.gro -r gros/nvt.gro -t bins/nvt.cpt -p topol.top -o bins/npt.tpr -po mdps/mdout.mdp") 
  system("gmx mdrun -deffnm bins/npt")
  system("mv bins/npt.log logs/npt.log") 
  system("mv bins/npt.gro gros/npt.gro") 
  system("echo 'Pressure' | gmx energy -f bins/npt.edr -o xvgs/pressure.xvg") 
  plotXVG("xvgs/pressure.xvg")
  system("echo 'Density' | gmx energy -f bins/npt.edr -o xvgs/density.xvg") 
  plotXVG("xvgs/density.xvg")
  
  ####### now for the production verion. Remove constraints => get rid of -r option
  ######## Hold Room Temp and Pressure and let the protein change shape ###########
  system("gmx grompp -f ../inputs/md500.mdp -c gros/npt.gro -t bins/npt.cpt -p topol.top -o bins/md.tpr -po mdps/mdout.mdp") 
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
  system("rm temp.top*")

  ################### Analyses ##############
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
}
setwd(base)
# NGLVieweR("EXT/pdbs/EXT.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("TRU/pdbs/TRU.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("BCC2_350/pdbs/BCC2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("BCC10_350/pdbs/BCC10.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("DIA2_350/pdbs/DIA2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
# NGLVieweR("FCC5_350/pdbs/FCC5.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("BCC2_500/pdbs/BCC2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("BCC10_500/pdbs/BCC10.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("DIA2_500/pdbs/DIA2.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")
NGLVieweR("FCC5_500/pdbs/FCC5.pdb")%>%addRepresentation("cartoon")%>%addRepresentation("ball+stick")

