from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open pdbs/CUB/LatCUB.pdb')
run(session, 'open pdbs/1L2Y_CA.pdb')
run(session, 'style #1:1,1912,1913,1929 sphere')
run(session, 'size #1:1912,1913 atomRadius 0.75')
run(session, 'size #1:1929 atomRadius 0.5')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:274')
run(session, 'align #2:1-3 to #1:1912,1913,1929 move chains')

run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-20 suppressBackboneDisplay false')
run(session, 'hide #2 target p')

run(session, 'view #2:1-20')  # who all is in the lattice slice zoomed in on
# 
# c = session.models[0]
# print(c.name, c.num_atoms, 'atoms', c.num_residues, 'residues', c.num_chains, 'chains')
# a = c.atoms[1928]
# print('Atom', a, 'at', a.coord)
# 
# t = session.models[2]
# print(t.name, t.num_atoms, 'atoms', t.num_residues, 'residues', t.num_chains, 'chains')

# 1912	1913
# 1912	1928
# 1912	2168
# 1913	1914
# 1913	1929
# 1913	2169
