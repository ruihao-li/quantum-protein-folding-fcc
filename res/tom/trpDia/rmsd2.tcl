color Element F green 
mol new pdbs/dia5.pdb
mol modstyle 0 0 VDW 0.300000 12.000000
mol modcolor 0 0 element
mol addrep 0
set diaStr "resid  500 502"; # makes sense since 1000 atoms total (5x5x5x8)
mol modselect 1 0 $diaStr
mol modstyle 1 0 VDW 0.900000 12.000000
mol modcolor 1 0 element
mol modmaterial 1 0 Transparent
set dia  [atomselect 0 "resid 1 to 1000"]
set d [$dia get {x y z}]
set dia2 [atomselect 0 $diaStr]


mol new pdbs/trpCen.pdb
mol modstyle 0 1 Tube 0.500000 12.000000
mol modcolor 0 1 ColorID 27
mol modmaterial 0 1 Transparent
mol addrep 1
mol modstyle 1 1 VDW 0.500000 12.000000
mol modcolor 1 1 ColorID 27
mol modmaterial 1 1 Transparent
set trp [atomselect 1 "all"]

translate to 0 0 -2
set file [open "outs/rmsd2.dat" w] 
for {set ip 1} {$ip < 20} {incr ip} {  
    set str "resid $ip to [expr $ip+1]"
    set trp2 [atomselect 1 $str] 
    set transformation_mat [measure fit $trp2 $dia2 weight {1 1}]
    $trp move $transformation_mat
    set t [$trp get {x y z}]
    puts $str
    set X {}
    for {set i 0} {$i < 20} {incr i} { 
	   set x {}
	   for {set j 0} {$j < 1000} {incr j} { 
	      lappend x [expr [veclength2 [vecsub [lindex $t $i] [lindex $d $j]]]]
	       }
	  set mn [tcl::mathfunc::min {*}$x]
	  puts $mn
	  lappend X $mn
	  display update
  }
  puts $X
  set mean [expr {[tcl::mathop::+ {*}$X 0.0] / max(1, [llength $X])}]
  set rmsd [expr sqrt($mean)]
  puts $file [list $ip $rmsd]
}
close $file

