color Element F green 
mol new pdbs/dia6.pdb
mol modstyle 0 0 VDW 0.300000 12.000000
mol modcolor 0 0 element
mol addrep 0
set diaStr "resid  693 1033"; 
mol modselect 1 0 $diaStr
mol modstyle 1 0 VDW 1.200000 12.000000
mol modcolor 1 0 element
mol modmaterial 1 0 Transparent

set dia2  [atomselect 0 $diaStr]
set diaA  [atomselect 0 "element N"]
set dA [$diaA get {x y z}]
set diaB  [atomselect 0 "element F"];# green atoms = B lattice
set dB [$diaB get {x y z}]

axes location Off;# use fiducials below instead
graphics top color red 
graphics top sphere {2.25 0 0} radius 0.1 resolution 10
graphics top color green 
graphics top sphere {0 2.25 0} radius 0.1 resolution 10
graphics top color blue 
graphics top sphere {0 0 2.25} radius 0.1 resolution 10

mol new pdbs/trpCen.pdb
mol modstyle 0 1 Tube 0.500000 12.000000
mol modcolor 0 1 ColorID 27
mol modmaterial 0 1 Transparent
mol addrep 1
mol modstyle 1 1 VDW 0.500000 12.000000
mol modcolor 1 1 ColorID 27
mol modmaterial 1 1 Transparent
set trp [atomselect 1 "all"]

set transRot1 [transabout {1 1 1} 1]

translate to 0 0 -4
set file [open "outs/rmsd6_AB.dat" w] 
# ip = 1 starts on F/green/B  and switches back and forth with each increase in ip   
for {set ip 1} {$ip < 20} {incr ip} {  
  set str "resid $ip to [expr $ip+1]"
  puts $str
  if {$ip%2 == 1} {
	  set ipd 1; # ip is odd (ipd is 1) when AA 1 is on green
	  } else {
	  set ipd 0
     }
  set trp2 [atomselect 1 $str] 
  set transformation_mat [measure fit $trp2 $dia2 weight {1 1}]
  $trp move $transformation_mat
  for {set theta 0} {$theta < 360} {incr theta} {  
    set t [$trp get {x y z}]
    set X {}
    # here i starts at 0 because lindex starts counting at zero 
    for {set i 0} {$i < 20} {incr i} { 
      set x {}
      if {$i%2 == $ipd} {
	  # so when ip=1 (and thus ipd=1) and i=0, we hit the else below (start on B/green)
      for {set j 0} {$j < 864} {incr j} { 
        lappend x [expr [veclength2 [vecsub [lindex $t $i] [lindex $dA $j]]]]
      } } else {
      for {set j 0} {$j < 864} {incr j} { 
        lappend x [expr [veclength2 [vecsub [lindex $t $i] [lindex $dB $j]]]]
	  }
     }
      set mn [tcl::mathfunc::min {*}$x]
      lappend X $mn
      display update
    }
    set mean [expr {[tcl::mathop::+ {*}$X 0.0] / max(1, [llength $X])}]
    set rmsd [expr sqrt($mean)]
    puts $file [list $ip $theta $rmsd]
    $trp move $transRot1
  }
}


close $file







