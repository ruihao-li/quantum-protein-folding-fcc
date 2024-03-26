from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open pdbs/LatCUB.pdb')
run(session, 'open pdbs/ethane.pdb')
# run(session, 'style #1:1,12,22,40,57,58,59, 61, 63,  sphere')
# run(session, 'size #1:57,58 atomRadius 0.75')
# run(session, 'size #1:12,22,40,59,61,63 atomRadius 0.5')
# run(session, 'size #1:1 atomRadius 0.5')
# run(session, 'view zalign #1:1 inFrontOf #1:2')
# run(session, 'align #2:1-8 to #1:57,58, 12,22,40, 59,61,63 move chains')

# 
run(session, 'style #1:2,10,12, 11,14, 13, 23, 15 sphere')
run(session, 'size  #1:2,10,12, 11,14, 13, 23, 15 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:2')
# run(session, 'align #2:1-8 to #1:11,14, 2,10,12,13, 23, 15  move chains')
run(session, 'align #2:1,2,3,6 to #1:11,14, 2,23  move chains')

run(session, 'open pdbs/cub.pdb')
# 
# run(session, 'show #2 target a')
# run(session, 'style #2 sphere')
# run(session, 'size #2 atomRadius 0.5')
# # # run(session, 'cartoon #2#1-8 suppressBackboneDisplay false')
# run(session, 'hide #2 target p')
# run(session, 'view #2:1-8')  # who all is in the lattice slice zoomed in on
