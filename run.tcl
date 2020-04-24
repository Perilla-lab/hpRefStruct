########################################################
# usage: proref ?-arg var? ...
#        -p pdb file
#        -m map file (in .map or .mrc format)
#        -r resolution of the map
#        -j job name
#        -n number of refinement iterations
#        -w weight of the density
#######################################################

source hprefstruct-2.tcl

set pdb #your pdb
set map #cryo-EM map

hprs -w 5.0 -p $pdb -m $map -j refine-1
