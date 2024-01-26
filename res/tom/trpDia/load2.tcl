color Element F green 
mol new pdbs/dia5.pdb

mol modstyle 0 0 VDW 0.300000 12.000000
mol modcolor 0 0 element

mol addrep 0
# 62*8 +1 = 497 = origin atom 1 of middle cube in 5x5x5, filling z first, then y, then x
# in first/unit cube, 4 is z=0 face, 6 is inner down and left of it, so 1=497 yields
set diaStr "resid  500 502"; # makes sense since 1000 atoms total (5x5x5x8)
mol modselect 1 0 $diaStr
mol modstyle 1 0 VDW 0.900000 12.000000
mol modcolor 1 0 element
mol modmaterial 1 0 Transparent


axes location Off;# use fiducials below instead
graphics top color red 
graphics top sphere {2.25 0 0} radius 0.1 resolution 10
graphics top color green 
graphics top sphere {0 2.25 0} radius 0.1 resolution 10
graphics top color blue 
graphics top sphere {0 0 2.25} radius 0.1 resolution 10


mol new pdbs/trpCen.pdb
set i 1
set str "resid $i to [expr $i+1]"
mol modstyle 0 1 Tube 0.500000 12.000000
mol modcolor 0 1 ColorID 27
mol modmaterial 0 1 Transparent

mol addrep 1
#mol modselect 1 1 $str
mol modstyle 1 1 VDW 0.500000 12.000000
mol modcolor 1 1 ColorID 27
mol modmaterial 1 1 Transparent

translate to 0 0 -2
set dia2 [atomselect 0 $diaStr]
set trp2 [atomselect 1 $str] 
set trp [atomselect 1 "all"]
set transformation_mat [measure fit $trp2 $dia2 weight {1 1}]
$trp move $transformation_mat

set t [$trp get {x y z}]
set dist2 [veclength2 [vecsub [lindex $t [expr $i-1]] [lindex $t $i]] ]
set dist [expr sqrt($dist2)]
puts "i is now $i and the distance is $dist"

