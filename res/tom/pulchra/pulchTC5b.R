#  https://github.com/euplotes/pulchra
# system("cc -O3 -o pulchra pulchra.c pulchra_data.c -lm")
# system("./pulchra")
# PULCHRA Protein Chain Restoration Algorithm version 3.06
# Usage: ./pulchra [options] <pdb_file>
#   The program default input is a PDB file.
# Output file <pdb_file.rebuild.pdb> will be created as a result.
# Valid options are:
#   
#   -v : verbose output (default: off)
# -n : center chain (default: off)
# -x : time-seed random number generator (default: off)
# -g : use PDBSG as an input format (CA=C-alpha, SC or CM=side chain c.m.)
# 
# -c : skip C-alpha positions optimization (default: on)
# -p : detect cis-prolins (default: off)
# -r : start from a random chain (default: off)
# -i pdbfile : read the initial C-alpha coordinates from a PDB file
# -t : save chain optimization trajectory to file <pdb_file.pdb.trajectory>
#   -u value : maximum shift from the restraint coordinates (default: 0.5A)
# 
# -e : rearrange backbone atoms (C, O are output after side chain) (default: off)
# -f : preserve initial coordinates (default: off, implies '-c' on and '-n' off)
# -b : skip backbone reconstruction (default: on)
# -q : optimize backbone hydrogen bonds pattern (default: off)
# -h : outputs hydrogen atoms (default: off)
# -s : skip side chains reconstruction (default: on)
# -o : don't attempt to fix excluded volume conflicts (default: on)
# -z : don't check amino acid chirality (default: on)
# system("./pulchra -v ../chiX/pdbs/FCC/fits/fcc2.pdb") # my best FCC = 1.363")
# to get pulchra to run from anywhere, you can soft link the binary into /usr/local/bin as follows
# cd /usr/local/bin
# sudo ln -s /Users/radivot/ccf/quantum/pulchra/pulchra pulchra


for (i in c("DIA","CUB","BCC","FCC")) 
  for (j in 2:19) {
     # i="DIA"
     # j=2
    system(paste0("pulchra -v ../TC5b/pdbs/",i,"/fits/",tolower(i),j,".pdb")) 
    system(paste0("mv ../TC5b/pdbs/",i,"/fits/",tolower(i),j,".rebuilt.pdb pdbTC5b/",i,j,".pdb")) 
  } 
