set x [lindex $argv 0]
foreach i [lindex $argv 1] {
set id [mol new target$x-$i.pdb]

[atomselect $id backbone] writepdb target$x-$i.bb.pdb
}

exit
