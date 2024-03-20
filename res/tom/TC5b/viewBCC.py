from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open pdbs/BCC/LatBCC.pdb')
run(session, 'open pdbs/1L2Y_CA.pdb')
run(session, 'style #1:1,2197,2198,2199 sphere')
run(session, 'size #1   atomRadius 0.25')
run(session, 'size #1:2197,2198,2199 atomRadius 0.75')
run(session, 'size #1:1 atomRadius 0.5')
run(session, 'view zalign #1:1 inFrontOf #1:533')
run(session, 'align #2:1-3 to #1:2197,2198,2199 move chains')
run(session, 'view #1:2197')
# RMSD between 3 atom pairs is 0.486 angstroms

run(session, 'style #1:1,2197,2198,2199,1860 sphere')
run(session, 'size #1:2197,2198,2199,1860 atomRadius 0.75')
run(session, 'align #2:1-4 to #1:2197,2198,2199,1860 move chains')
run(session, 'view #1:2197,2198,2199,1860')
# RMSD between 4 atom pairs is 0.766 angstroms

run(session, 'style #1:1,2197,2198,2199,1860,2225 sphere')
run(session, 'size #1:2197,2198,2199,1860,2225 atomRadius 0.75')
run(session, 'align #2:1-5 to #1:2197,2198,2199,1860,2225 move chains')
run(session, 'view #1:2197,2198,2199,1860,2225') 
# RMSD between 5 atom pairs is 0.916 angstroms


run(session, 'style #1:1,2197,2198,2199,1860,2225,1862 sphere')
run(session, 'size #1:2197,2198,2199,1860,2225,1862 atomRadius 0.75')
run(session, 'align #2:1-6 to #1:2197,2198,2199,1860,2225,1862 move chains')
run(session, 'view #1:2197,2198,2199,1860,2225,1862') 
# RMSD between 6 atom pairs is 1.257 angstroms

run(session, 'style #1:1,2197,2198,2199,1860,2225,1862,1889 sphere')
run(session, 'size #1:2197,2198,2199,1860,2225,1862,1889 atomRadius 0.75')
run(session, 'align #2:1-7 to #1:2197,2198,2199,1860,2225,1862,1889 move chains')
run(session, 'view #1:2197,2198,2199,1860,2225,1862,1889') 
# RMSD between 7 atom pairs is 1.291 angstroms


run(session, 'show #2 target a')
run(session, 'style #2 sphere')
run(session, 'size #2 atomRadius 0.5')
run(session, 'cartoon #2#1-20 suppressBackboneDisplay false')
run(session, 'hide #2 target p')


c = session.models[0]
print(c.name, c.num_atoms, 'atoms', c.num_residues, 'residues', c.num_chains, 'chains')
a = c.atoms[2196]
print('Atom', a, 'at', a.coord)
# 
# 2197	2198
# 2198	2199
# 2198	2223
# 2198	2225
# 2198	2535
# 2198	2537
# 2198	2561
# 2198	2563
