mol new pdbs/diaR2p3C.pdb
mol modstyle 0 0 VDW 0.300000 12.000000
mol modcolor 0 0 element

#3x3x3 cube is initially 12x12x12. Centering cyl at (5,6,0) => 2.3*5 = 11.5 and 2.3*6=13.8, 2.3*12=27.6
graphics top color red 
graphics top cylinder {11.5 13.8 0} {11.5 13.8 27.6} radius 2.3 resolution 60 filled no 
set cylCen {11.5 13.8 13.8};# center of 4-barrel cyl

#axes location origin
axes location Off;# use fiducials below instead
graphics top sphere {2.3 0 0} radius 0.1 resolution 10
graphics top color green 
graphics top sphere {0 2.3 0} radius 0.1 resolution 10
graphics top color blue 
graphics top sphere {0 0 2.3} radius 0.1 resolution 10

####### learn from logfile console






