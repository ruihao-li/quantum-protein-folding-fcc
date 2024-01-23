set transRot1 [transabout {0 0 1} 1]

for {set i 0} {$i < 360} {incr i} {  
$alf move $transRot1
$trp move $transRot1
display update
}

#set transRot [trans axis z 90 ]
#set transRot [trans axis z 45 ]

set transRot180 [transabout {0 0 1} 180]
if (0) {
$trp move $transRot180
}

#set transRot30 [transabout {0 0 1} 30]
#$trp move $transRot30
#display update
#sleep 2
#$trp move $transRot30
#display update

#puts $transShift
#graphics top sphere {2.3 0 0} radius 0.1 resolution 10
#set transShift [transvec {-5 -6 0}]
#puts $transShift
#$trp move $transShift

#set transRot [transaxis z 90]
#puts $transRot
#$trp move $transRot

#mol fix 0
#mol fix 1
#rotate z by 360 5
##mol free 0
#mol free 1





