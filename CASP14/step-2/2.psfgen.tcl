set ip 46
mol new refine-1_run4_0001.pdb
#target$ip.resaved.pdb
set all [atomselect top all]
$all set segname PA
[atomselect top "resname HIS"] set resname HSD
#[atomselect top "resname HSD and resid 36"] set resname HSP
#[atomselect top "resname HSD and resid 62"] set resname HSE
$all writepdb target$ip.tmp.pdb
mol delete all

package require psfgen

topology /misc/chaoyi/toppar/top_all36_prot.rtf
pdbalias atom ILE CD1 CD
#pdbalias residue HIS HSD
segment PA {
    first NTER
    last CTER
    pdb target$ip.tmp.pdb
}
coordpdb target$ip.tmp.pdb PA
guesscoord

writepsf target$ip.psf
writepdb target$ip.pdb

mol delete all
resetpsf

package require solvate
solvate target$ip.psf target$ip.pdb -o target$ip.solvated -s S -b 2.4 -t 10

package require autoionize
autoionize -psf target$ip.solvated.psf -pdb target$ip.solvated.pdb -sc 0.15 -seg ION -o target$ip.preped

mol delete all
resetpsf

set id [mol new target$ip.preped.psf]
mol addfile target$ip.preped.pdb
set all [atomselect $id all]
$all moveby [vecinvert [measure center $all weight mass]]
$all moveby [vecinvert [measure center $all weight mass]]
$all moveby [vecinvert [measure center $all weight mass]]
$all writepdb target$ip.preped.pdb

file delete target$ip.solvated.psf target$ip.solvated.pdb
#file delete target$ip.tmp.pdb target$ip.solvated.psf target$ip.solvated.pdb

exit
