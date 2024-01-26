incr i

set str "resid $i to [expr $i+1]"

mol modselect 1 0 $diaStr
mol modstyle 1 0 VDW 0.900000 12.000000
mol modcolor 1 0 element
mol modmaterial 1 0 Transparent

translate to 0 0 -2
set dia2 [atomselect 0 $diaStr]
set trp2 [atomselect 1 $str] 
set transformation_mat [measure fit $trp2 $dia2 weight {1 1}]
$trp move $transformation_mat

set t [$trp get {x y z}]
set dist2 [veclength2 [vecsub [lindex $t [expr $i-1]] [lindex $t $i]] ]
set dist [expr sqrt($dist2)]
puts "i is now $i and the distance is $dist"

