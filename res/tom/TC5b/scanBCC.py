# Bi-directional version of rmsd.py
import time
import re
import json
import numpy as np
from chimerax.core.commands import run  
run(session, 'close')
run(session, 'open pdbs/BCC/LatBCC.pdb')
run(session, 'open pdbs/1L2Y_CA.pdb')
run(session, 'size #1   atomRadius 0.25')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'style #1:1,2197,2198,2199 sphere')
run(session, 'size #1:2197,2198,2199 atomRadius 0.75')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:533')
run(session, 'align #2:1-3 to #1:2197,2198,2199 move chains')
# run(session, 'align #2:1-4 to #1:2197,2198,2199,1860 move chains')
# run(session, 'view #1:2197')
# run(session, 'align #2:1-4 to #1:2197,2198,2199,1860 move chains')



S=list()
f = open("pdbs/BCC/conBCC.txt", "r")
import re
for x in f:
  y=x.strip()
  S.append(re.split('\s+', y))
f.close()

L=[list() for i in range(0,4395)]

# print(S[0])
# print(S[1])
# print(S[1][1])
# print(S[1][0])

# print(L)
# print(L[1])


for i in range(1,len(S)):
   # print(i)
   L[int(S[i][0])].append(S[i][1])
   L[int(S[i][1])].append(S[i][0])


d={}
for i in range(1,4395):
  d[str(i)]=L[i]

R=[0 for J in range(2,20)]
B1=[0 for J in range(2,20)]
# #J in [2, 19] is midpoint of starting triplet. J first runs up to 19 and then down to 2
# J=19 # starting at the other end gives a worse fit, RMSD of 3.4
# J=10 # starting in the middle gives 3.05, not bad
# J=2  # 2.9 as in rmsd.py
# run(session, 'align #2:1-4 to #1:2197,2198,2199,1860 move chains')
for J in range(2,20):
  print("J is:",J)
  # base1=["689","690","691"]  #Diamond
  # base1=["1912","1913","1929"] # Cubic
  # base1=["2661","2662","2663"]  #FCC
  base1=["2197","2198","2199"]  #BCC
  # base1=["2198","2199","1860"]  #BCC
  base2=[str(J-1),str(J),str(J+1)]
  base1S=', '.join(base1)
  base2S=', '.join(base2)
  print(base2S)
  RMSD=[0 for i in range(0,18)]
  # RMSD[J-2]=np.around(run(session, 'align #2:'+base2S+' to #1:'+base1S+' cutoffDistance 3 move chains')[2],decimals=3)
  RMSD[J-2]=np.around(run(session, 'align #2:'+base2S+' to #1:'+base1S+' move chains')[2],decimals=3)
  print(RMSD)

  print(base1)
  print(base1[-1])
  for j in range(J+2,21):
    print(j)
    rmsd={}
    nxt=d[base1[-1]] # grab last one while growing to larger values
    print("\n\nj is:",j)
    base2.append(str(j))
    base2S=', '.join(base2)
    for i in nxt:
      if i not in base1:
        tmp1=base1.copy()
        tmp1.append(i)
        tmp1S=', '.join(tmp1)
        rmsd[i]=run(session, 'align #2:'+base2S+' to #1:'+tmp1S+' move chains')[2]
    kys=list(rmsd.keys())
    vals=list(rmsd.values())
    minpos = vals.index(min(vals))
    ky=kys[minpos]
    base1.append(ky)
    base1S=', '.join(base1)
    RMSD[j-3]=np.around(rmsd[ky],decimals=3)

  for j in range(J-2,0,-1):
    rmsd={}
    nxt=d[base1[0]] # grab first one now instead of last, to grow toward smaller AA #s
    print("\n\nj is:",j)
    base2.insert(0,str(j)) # first arg is index, so insert into front
    base2S=', '.join(base2)
    for i in nxt:
      if i not in base1:
        tmp1=base1.copy()
        tmp1.insert(0,i)
        tmp1S=', '.join(tmp1)
        rmsd[i]=run(session, 'align #2:'+base2S+' to #1:'+tmp1S+' move chains')[2]
    kys=list(rmsd.keys())
    vals=list(rmsd.values())
    minpos = vals.index(min(vals))
    ky=kys[minpos]
    base1.insert(0,ky)
    base1S=', '.join(base1)
    RMSD[j-1]=np.around(rmsd[ky],decimals=3)

  run(session, 'align #2:'+base2S+' to #1:'+base1S+' move chains')[2] #make final move
  print("\n\nRMSD vs #of AA starting at 4:\n",RMSD,"\n")
  R[J-2]=max(RMSD)
  B1[J-2]=base1S
  print("Ending J=",J," we have R=",R)

print(B1)
f = open("pdbs/BCC/BCCfits.txt", "w")
for x in B1:
  f.write(x+'\n')
f.close()


run(session, 'style #1:1,'+base1S+' sphere')
run(session, 'size #1:'+base1S+' atomRadius 0.75')
run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-20 suppressBackboneDisplay false')
run(session, 'hide #2 target p')
run(session, 'view #2:'+base2S)
print(R)

# BCC lattic fit                              X     X
# 1.755, 1.733, 1.919, 1.744, 1.586, 1.632, 1.53, 1.53, 1.746, 2.392, 1.746, 1.746, 1.656, 1.746, 1.848, 1.656, 1.651, 1.848] # 2JOF
# 8 (2^3) basis vecs  (density=0.68 is close to maximum)
#    X     3      4                            8     9     10     11     12                         X
# 1.441, 1.653, 1.851, 1.851, 1.717, 1.703, 1.703, 1.522, 1.868, 1.761, 1.75, 1.75, 1.695, 1.767, 2.108, 1.894, 1.666, 1.887] # 1L2Y


# FCC lattic fit          X                   X
# 1.363, 1.387, 1.432, 1.377, 1.797, 1.538, 1.38, 1.52, 1.418, 1.447, 1.513, 1.641, 1.507, 1.651, 1.449, 1.595, 1.565, 1.653] # 2JOF
# latfit gave  cRMSD = 1.4113 Angstroms
# 12 basis vecs    (density=0.76, maximum)
#                         X
# 1.646, 1.621, 1.439, 1.212, 1.239, 1.477, 1.489, 1.296, 1.585, 1.326, 1.714, 1.328, 1.583, 1.491, 1.6, 1.554, 1.575, 1.699] # 1L2Y

# cubic lattic fit                             X
# [1.844, 1.929, 1.871, 1.871, 1.871, 1.745, 1.713, 1.895, 2.845, 2.845, 2.845, 2.125, 2.025, 2.025, 2.025, 2.025, 2.025, 2.025] # 2JOF
# latfit gave  cRMSD = 1.9238 Angstroms
# 6 basis vecs  (density=0.52)
##                       X
# 2.273, 2.303, 2.312, 1.928, 2.312, 2.544, 1.982, 1.982, 2.508, 1.982, 1.982, 1.982, 1.982, 2.302, 2.227, 2.227, 2.648, 2.227] # 1L2Y


# diamond   X                   |<------------ same structure ------------------->|
# [2.893, 2.829, 2.989, 2.989, 3.047, 3.047, 3.047, 3.047, 3.047, 3.047, 3.047, 3.047, 3.593, 3.632, 3.632, 3.047, 3.437, 3.437] # 2JOF
# 4 (2^2) basis vecs (density=0.34)
#   X      X      X            X      X
#2.991, 2.991, 2.991, 3.671, 2.991, 2.991, 3.711, 3.711, 3.509, 3.827, 3.827, 3.827, 3.827, 3.827, 3.105, 3.105, 3.436, 3.436]  # 1L2Y





#
# savePDB=True
# savePDB=False
# if savePDB:
#   run(session, 'select #1:'+base1S)
#   run(session, 'save pdbs/trpCageDia.pdb format pdb models sel selectedOnly true')
#   run(session, 'select #2')
#   run(session, 'save pdbs/trpCage.pdb format pdb models sel selectedOnly true')
