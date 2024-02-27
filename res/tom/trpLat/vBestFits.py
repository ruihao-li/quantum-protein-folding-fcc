from chimerax.core.commands import run  
from chimerax.std_commands import align  
run(session, 'close')
run(session, 'open ../trpDia/pdbs/trpCen.pdb')
run(session, 'show #1 target a')
run(session, 'style #1 sphere')
run(session, 'size #1 atomRadius 0.5')
run(session, 'cartoon #1#1-20 suppressBackboneDisplay false')
run(session, 'hide #1 target p')
run(session, 'view #1')

run(session, 'open pdbs/FCC/fits/fcc2.pdb') # my best FCC = 1.363
run(session, 'style #2 sphere')
run(session, 'color #2 blue')
run(session, 'size #2 atomRadius 0.25')
run(session, 'align #2 to #1 move chains')[2]

run(session, 'open pdbs/CUB/fits/cub8.pdb') # my best CUB = 1.713
run(session, 'style #3 sphere')
run(session, 'color #3 orange')
run(session, 'size #3 atomRadius 0.25')
run(session, 'align #3 to #1 move chains')[2]

run(session, 'open pdbs/FCC/LatFitFCC.pdb') # latfit's FCC = 1.4113
run(session, 'delete #4/P')
run(session, 'style #4/L sphere')
run(session, 'color #4/L blue')
run(session, 'transparency #4/L 55 target c')
run(session, 'size #4/L atomRadius 0.25')
run(session, 'align #4/L to #1 move chains')[2]  # 1.345

run(session, 'open pdbs/CUB/LatFitCUB.pdb') # latfit's CUB = 1.9238
run(session, 'delete #5/P')
run(session, 'style #5/L sphere')
run(session, 'size #5/L atomRadius 0.25')
run(session, 'transparency #5/L 55 target c')
run(session, 'color #5/L orange')
run(session, 'align #5/L to #1 move chains')[2]  # 1.760 

run(session, 'open pdbs/BCC/fits/bcc7.pdb') # my best BCC = 1.53
run(session, 'style #6 sphere')
run(session, 'color #6 green')
run(session, 'size #6 atomRadius 0.25')
run(session, 'align #6 to #1 move chains')[2]

run(session, 'open pdbs/DIA/fits/dia3.pdb') # my best DIA = 2.829
run(session, 'style #7 sphere')
run(session, 'color #7 red')
run(session, 'size #7 atomRadius 0.25')
run(session, 'align #7 to #1 move chains')[2]
