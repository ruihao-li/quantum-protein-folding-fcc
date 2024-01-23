# this takes some time to run (guess ~15 minutes) 
set d216  [atomselect 0 "resid 1 to 216"]
set d [$d216 get {x y z}]
set t11   [atomselect 2 "resid 10 to 20"]

set file [open "outs/theta.dat" w] 
set transRot1 [transabout {0 0 1} 1]

for {set theta 0} {$theta < 360} {incr theta} {  
  set t [$t11 get {x y z}]
  set X {}
  for {set i 0} {$i < 11} {incr i} { 
    set x {}
    for {set j 0} {$j < 216} {incr j} { 
     lappend x [expr [veclength2 [vecsub [lindex $t $i] [lindex $d $j]]]]
    }
    set min [tcl::mathfunc::min {*}$x]
    lappend X $min
  }
  puts [list $theta]
  set mean [expr {[tcl::mathop::+ {*}$X 0.0] / max(1, [llength $X])}]
  set rmsd [expr sqrt($mean)]
  puts $file [list $theta $rmsd]
  $t11 move $transRot1
#  $trp move $transRot1
#  display update
}

close $file

