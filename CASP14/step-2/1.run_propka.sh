
pdb2pqr --ff=charmm --apbs-input --ph-calc-method=propka --with-ph=7.4 -v --summary  $1 $1.pqr

grep "HIS "  $1.pqr > $1.his

echo "output is  $1.pqr"
