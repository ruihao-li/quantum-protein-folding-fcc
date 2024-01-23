set dia [atomselect 0 "all"]
set alf [atomselect 1 "all"]
set alf9 [atomselect 1 "resid 6 to 14"]
set trp9 [atomselect 2 "resid 1 to 9"]
set trp [atomselect 2 "all"]

set transformation_mat [measure fit $trp9 $alf9]
$trp move $transformation_mat

