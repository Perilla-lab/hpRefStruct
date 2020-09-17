x=$1
sed "s/set x baba/set x $x/" mini.namd > mini.namd2
/home/chaoyi/local/bin/namd2 +p4 mini.namd2 > 1.mini.log &
