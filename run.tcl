########################################################
# usage: proref ?-arg var? ...
#        -p pdb file
#        -m map file (in .map or .mrc format)
#        -r resolution of the map
#        -j job name
#        -n number of refinement iterations
#        -w weight of the density
#######################################################

source hprefstruct.tcl

set pdb [glob ????.pdb]
set map [glob emd_*.map]

proref -w 5.0 -p $pdb -m $map -j ref-1
