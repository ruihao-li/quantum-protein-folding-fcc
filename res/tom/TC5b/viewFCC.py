from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open pdbs/FCC/LatFCC.pdb')
run(session, 'open pdbs/1L2Y_CA.pdb')
run(session, 'style #1:1,2661,2662,2663 sphere')
run(session, 'size #1   atomRadius 0.25')
run(session, 'size #1:2661,2662 atomRadius 0.75')
run(session, 'size #1:2663 atomRadius 0.5')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:533')
run(session, 'align #2:1-3 to #1:2661,2662,2663 move chains')

run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-20 suppressBackboneDisplay false')
run(session, 'hide #2 target p')

run(session, 'view #2:1-20')

# c = session.models[0]
# print(c.name, c.num_atoms, 'atoms', c.num_residues, 'residues', c.num_chains, 'chains')
# a = c.atoms[2660]
# print('Atom', a, 'at', a.coord)


# 2661	2662
# 2661	2663
# 2661	2664
# 2662	2663
# 2662	2664
# 2662	2665
# 2662	2667
# 2662	2705
# 2662	2708
# 2662	2709
