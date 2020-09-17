file mkdir results
set x [lindex $argv 0]
set i 1
cd 1
    set id [mol new target$x.preped.psf]
	mol addfile output/mini-1.coor
	[atomselect $id all] set occupancy 1.0
	[atomselect $id "resname HSE HSD HSP"] set resname HIS
	[atomselect $id "resname ILE and name CD"] set name CD1
	[atomselect $id "protein and noh"] writepdb ../results/target$x-${i}.pdb

exit
