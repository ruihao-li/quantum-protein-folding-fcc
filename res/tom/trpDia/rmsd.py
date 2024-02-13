import time
import re
import json
import numpy as np
from chimerax.core.commands import run  
run(session, 'close')
run(session, 'open pdbs/LatDia6.pdb')
run(session, 'open pdbs/trpCen.pdb')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:2')
base1=["689","690","691"]
base2=["1","2","3"]
base1S=', '.join(base1)
base2S=', '.join(base2)
RMSD=[0 for i in range(0,18)]
RMSD[0]=np.around(run(session, 'align #2:'+base2S+' to #1:'+base1S+' move chains')[2],decimals=3)
print(RMSD)
S=list()
f = open("pdbs/con.txt", "r")
import re
for x in f:
  y=x.strip()
  S.append(re.split('\s+', y))
f.close()

L=[list() for i in range(0,1729)]

for i in range(1,len(S)):
   L[int(S[i][0])].append(S[i][1])
   L[int(S[i][1])].append(S[i][0])

d={}
for i in range(1,1729):
  d[str(i)]=L[i]
print(base1)
print(base1[-1])
for j in range(4,21):
  rmsd={}
  nxt=d[base1[-1]]
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

run(session, 'align #2:'+base2S+' to #1:'+base1S+' move chains')[2] #make final move
print("\n\nRMSD vs #of AA starting at 4:\n",RMSD,"\n")

run(session, 'style #1:1,'+base1S+' sphere')
run(session, 'size #1:'+base1S+' atomRadius 0.75')
run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-20 suppressBackboneDisplay false')
run(session, 'hide #2 target p')
run(session, 'view #2:'+base2S)

savePDB=False
savePDB=True
if savePDB:
  run(session, 'select #1:'+base1S)
  run(session, 'save pdbs/trpCageDia.pdb format pdb models sel selectedOnly true')
  run(session, 'select #2')
  run(session, 'save pdbs/trpCage.pdb format pdb models sel selectedOnly true')
