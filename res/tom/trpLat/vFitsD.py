from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open ../trpDia/pdbs/trpCen.pdb')
run(session, 'show #1 target a')
run(session, 'style #1 sphere')
run(session, 'size #1 atomRadius 0.5')
run(session, 'cartoon #1#1-20 suppressBackboneDisplay false')
run(session, 'hide #1 target p')
run(session, 'view #1')

for i in range(0,18):
  run(session, 'open pdbs/DIA/fits/dia'+str(i+2)+'.pdb')
  run(session, 'style #'+str(i+2)+' sphere')
  run(session, 'size #'+str(i+2)+' atomRadius 0.25')
  run(session, 'align #'+str(i+2)+' to #1 move chains')[2] #make final move

# ## chunck proves that dia6 to dia13 are all the same
# for i in range(0,8):
#   run(session, 'open pdbs/DIA/fits/dia'+str(i+6)+'.pdb')
#   run(session, 'style #'+str(i+2)+' sphere')
#   run(session, 'size #'+str(i+2)+' atomRadius 0.25')
#   run(session, 'align #'+str(i+2)+' to #1 move chains')[2] #make final move
