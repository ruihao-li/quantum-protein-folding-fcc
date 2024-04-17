from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open pdbs/5AWL_CA.pdb')
run(session, 'show #1 target a')
run(session, 'style #1 sphere')
run(session, 'size #1 atomRadius 0.5')
run(session, 'cartoon #1#1-10 suppressBackboneDisplay false')
run(session, 'hide #1 target p')
run(session, 'view #1')

run(session, 'open pdbs/DIA/dia2.pdb') 
run(session, 'style #2 sphere')
run(session, 'color #2 red')
run(session, 'size #2 atomRadius 0.25')
run(session, 'align #2 to #1 move chains')[2]

run(session, 'open pdbs/FCC/fcc4.pdb') 
run(session, 'style #3 sphere')
run(session, 'color #3 blue')
run(session, 'size #3 atomRadius 0.25')
run(session, 'align #3 to #1 move chains')[2]
# run(session, 'cartoon suppressBackboneDisplay false')



