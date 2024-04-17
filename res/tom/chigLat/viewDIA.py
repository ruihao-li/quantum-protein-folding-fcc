from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open ../TC5b/pdbs/DIA/LatDIA.pdb')
run(session, 'open pdbs/5awl_CA.pdb')
# run(session, 'style #1:1,1033,1034,1035,1030 sphere')
# run(session, 'size #1:1033,1034 atomRadius 0.75')
# run(session, 'size #1:1035,1030 atomRadius 0.5')
# run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:2')
run(session, 'align #2:1-4 to #1:1033,1034,1035,1030 move chains')

run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-10 suppressBackboneDisplay false')
run(session, 'hide #2 target p')
run(session, 'view #2:1-10')

# c = session.models[0]
# print(c.name, c.num_atoms, 'atoms', c.num_residues, 'residues', c.num_chains, 'chains')
# a = c.atoms[1032]
# print('Atom', a, 'at', a.coord)
# a = c.atoms[1033]
# print('Atom', a, 'at', a.coord)
# 
# t = session.models[2]
# print(t.name, t.num_atoms, 'atoms', t.num_residues, 'residues', t.num_chains, 'chains')
