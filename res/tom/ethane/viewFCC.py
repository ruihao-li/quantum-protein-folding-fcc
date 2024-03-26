from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open pdbs/LatFCC.pdb')
run(session, 'open pdbs/ethane.pdb')
run(session, 'style #1:1 sphere')
run(session, 'size #1 atomRadius 0.1')
run(session, 'size #1:7,29, 6,28 atomRadius 0.5')
run(session, 'view zalign #1:7 inFrontOf #1:29')
# run(session, 'align #2:1,2,3,6 to #1:7,29, 6,28  move chains')
run(session, 'align #2:1,2,3,4,5,6,7,8 to #1:7,29,6,21,2,28,32,15  move chains')
# 
# run(session, 'show #2 target a')
# run(session, 'style #2 sphere')
# run(session, 'size #2 atomRadius 0.5')
# run(session, 'cartoon #2#1-20 suppressBackboneDisplay false')
# run(session, 'hide #2 target p')
# 
# run(session, 'view #2:1-8')

run(session, 'open pdbs/fcc.pdb')
