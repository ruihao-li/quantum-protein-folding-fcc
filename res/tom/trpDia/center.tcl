graphics 0 delete all

set transShift [transoffset {-11.5 -13.8 0}]
$trp move $transShift
$dia move $transShift
$alf move $transShift

translate by 1.1 1.1 0.3;# these are screen coords!!!
mouse mode center;#=typing c.  click on big Sulfur (yellow) 2-3 times, then r to start rotating about center
