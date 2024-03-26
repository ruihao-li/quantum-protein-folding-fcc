from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open pdbs/LatDIA.pdb')
run(session, 'open pdbs/ethane.pdb')
run(session, 'style #1:1,12,22,40,57,58,59, 61, 63,  sphere')
run(session, 'size #1:57,58 atomRadius 0.75')
run(session, 'size #1:12,22,40,59,61,63 atomRadius 0.5')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:2')
run(session, 'align #2:1-8 to #1:57,58, 12,22,40, 59,61,63 move chains')
# run(session, 'align #2:1,2,3 to #1:57,58, 59 move chains')
# run(session, 'align #2:1-4 to #1:40, 57,58, 59 move chains')

run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
# run(session, 'hide #2 target p')
# run(session, 'view #2:1-8')

run(session, 'open pdbs/dia.pdb')
