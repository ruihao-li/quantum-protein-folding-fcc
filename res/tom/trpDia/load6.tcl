color Element F green 
mol new pdbs/dia6.pdb

mol modstyle 0 0 VDW 0.300000 12.000000
mol modcolor 0 0 element

mol addrep 0
# unit cell is
#~ ATOM      1  C   DIA B   1       0.000   0.000   0.000  1.00  0.00           N  
#~ ATOM      2  C   DIA B   2       0.000   2.000   2.000  1.00  0.00           N  
#~ ATOM      3  C   DIA B   3       2.000   0.000   2.000  1.00  0.00           N  
#~ ATOM      4  C   DIA B   4       2.000   2.000   0.000  1.00  0.00           N  
#~ ATOM      5  C   DIA B   5       3.000   3.000   3.000  1.00  0.00           F  
#~ ATOM      6  C   DIA B   6       3.000   1.000   1.000  1.00  0.00           F  
#~ ATOM      7  C   DIA B   7       1.000   3.000   1.000  1.00  0.00           F  
#~ ATOM      8  C   DIA B   8       1.000   1.000   3.000  1.00  0.00           F


# (2*36+2*6+2)*8 +5 = 693    6x6x6 made filling z first, then y, then x
# (3*36+3*6+3)*8 +1 = 1033 
set diaStr "resid  693 1033"; 
mol modselect 1 0 $diaStr
mol modstyle 1 0 VDW 1.200000 12.000000
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
