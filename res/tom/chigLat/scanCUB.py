# Bi-directional version of rmsd.py
import time
import re
import json
import numpy as np
from chimerax.core.commands import run  
run(session, 'close')
run(session, 'open ../TC5b/pdbs/CUB/LatCUB.pdb')
run(session, 'open pdbs/5AWL_CA.pdb')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:274')

run(session, 'style #1:1,1912,1913,1929 sphere')
run(session, 'size #1:1912,1913,1929 atomRadius 0.75')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:274')
run(session, 'align #2:1-3 to #1:1912,1913,1929 move chains')



S=list()
f = open("../TC5b/pdbs/CUB/conCUB.txt", "r")
import re
for x in f:
  y=x.strip()
  S.append(re.split('\s+', y))
f.close()

L=[list() for i in range(0,4097)]

print(S[0])
print(S[1])
print(S[1][1])
print(S[1][0])

print(L)
print(L[1])
 

for i in range(1,len(S)):
   print(i)
   L[int(S[i][0])].append(S[i][1])
   L[int(S[i][1])].append(S[i][0])


d={}
for i in range(1,4097):
  d[str(i)]=L[i]

R=[0 for J in range(2,10)]
B1=[0 for J in range(2,10)]
for J in range(2,10):
  print("J is:",J)
  base1=["1912","1913","1929"]
  base2=[str(J-1),str(J),str(J+1)]
  base1S=', '.join(base1)
  base2S=', '.join(base2)
  RMSD=[0 for i in range(0,8)]
  RMSD[J-2]=np.around(run(session, 'align #2:'+base2S+' to #1:'+base1S+' move chains')[2],decimals=3)
  print(RMSD)
  
  print(base1)
  print(base1[-1])
  for j in range(J+2,11):
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
f = open("pdbs/CUB/CUBfits.txt", "w")
for x in B1:
  f.write(x+'\n')
f.close()


run(session, 'style #1:1,'+base1S+' sphere')
run(session, 'size #1:'+base1S+' atomRadius 0.75')
run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-10 suppressBackboneDisplay false')
run(session, 'hide #2 target p')
run(session, 'view #2:'+base2S)
print(R)
