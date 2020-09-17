file mkdir results
set x 46
foreach i " 1 2 3 4 5 6 7 8 9 10 11 12" {
    for {set j 0} {$j < 10} {incr j} {
	set id [mol new ../target${x}.preped.psf]
	mol addfile run-${i}/SA-${j}.coor
	[atomselect $id all] set occupancy 1.0
	[atomselect $id "resname HSE HSD HSP"] set resname HIS
	[atomselect $id "resname ILE and name CD"] set name CD1
	[atomselect $id protein] writepdb results/target${x}-${i}${j}.pdb
	mol delete $id
    }
}
exit
