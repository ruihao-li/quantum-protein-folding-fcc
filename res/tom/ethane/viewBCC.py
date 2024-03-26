from chimerax.core.commands import run
from chimerax.std_commands import align
run(session, 'close')
run(session, 'open pdbs/LatBCC.pdb')
run(session, 'open pdbs/ethane.pdb')
run(session, 'style #1:3,5,9,2,15,8,12,14  sphere')
run(session, 'size #1:3,5,9,2,15,8,12,14 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:2')
# run(session, 'align #2:1-8 to #1:3,5,9,2,15,8,12,14 move chains')
# run(session, 'align #2:1-8 to #1:2,15,3,5,9,8,12,14 move chains')
run(session, 'align #2:1,2,3,6 to #1: 2,15,9, 8  move chains')

run(session, 'open pdbs/bcc.pdb')

# run(session, 'show #2 target a')
# run(session, 'style #2 sphere')
# run(session, 'size #2 atomRadius 0.5')

