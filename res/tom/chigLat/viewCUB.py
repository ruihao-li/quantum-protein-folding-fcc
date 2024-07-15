from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open ../TC5b/pdbs/CUB/LatCUB.pdb')
run(session, 'open pdbs/5awl_CA.pdb')
run(session, 'style #1:1,1912,1913,1929 sphere')
run(session, 'size #1:1912,1913 atomRadius 0.75')
run(session, 'size #1:1929 atomRadius 0.5')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:274')
run(session, 'align #2:1-3 to #1:1912,1913,1929 move chains')

run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-10 suppressBackboneDisplay false')
run(session, 'hide #2 target p')

run(session, 'view #2:1-10')  
