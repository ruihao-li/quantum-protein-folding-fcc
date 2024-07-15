import time
import re
import json
import numpy as np
from chimerax.core.commands import run  
run(session, 'close')
# run(session, 'open pdbChig/EXT.pdb')   # 0.2  30
run(session, 'open pdbChig/TRU.pdb')   # 2.8  5.5
# run(session, 'open pdbChig/FCC4.pdb')  # 1.3  6.8
# run(session, 'open pdbChig/DIA2.pdb')  # 1.33 7.5
run(session, 'delete solvent')
#       DIA2
# DIA   1.88,  1.936, 2.529, 2.273, 2.273, 2.135, 2.135, 2.135  RMSD
# DIA   7.462, 7.462, 3.897, 7.462, 7.462, 7.462, 7.462, 7.462  end-to-end
# DIA   1.331, 1.322, 0.834, 0.628, 0.603, 0.78,  0.78,  0.851  # of Hbonds

#                     FCC4
# FCC   1.421, 1.431, 1.104, 1.35,  1.424, 1.144, 1.582, 1.596  RMSD
# FCC   7.6,   3.8,   6.582, 6.582, 6.582, 5.374, 3.8,   7.6    end-to-end
# FCC   2.025, 2.054, 1.279, 2.432, 1.98,  1.716, 1.182, 1.194  # of Hbonds


def H(r):
  return (1-pow(r/4,6))/(1-pow(r/4,8))

Htot=[0 for J in range(2,10)]
Dee=[0 for J in range(2,10)]

r1=run(session, 'distance #1:1@o #1:10@n')
print(r1)

r2=run(session, 'distance #1:3@n #1:8@o')
print(r2)

r3=run(session, 'distance #1:3@o #1:7@n')
print(r3)

h=H(r1)+H(r2)+H(r3)

x=run(session, 'distance #1:1@ca #1:10@ca')
print(x)
print(type(x))
Dee[0]=np.around(x,decimals=3)
Htot[0]=np.around(h,decimals=3)
print(Dee)
print(Htot)

lat="BCC"
lat="DIA"
lat="FCC"

for J in range(2,10):
  print("J is:",J)
  run(session, 'close')
  run(session, 'open pdbChig/'+lat+str(J)+'.pdb') 
  r1=run(session, 'distance #1:1@o #1:10@n')
  print(r1)
  r2=run(session, 'distance #1:3@n #1:8@o')
  print(r2)
  r3=run(session, 'distance #1:3@o #1:7@n')
  print(r3)
  h=H(r1)+H(r2)+H(r3)
  x=run(session, 'distance #1:1@ca #1:10@ca')
  print(x)
  Dee[J-2]=np.around(x,decimals=3)
  Htot[J-2]=np.around(h,decimals=3)
print(Dee)
print(Htot)








