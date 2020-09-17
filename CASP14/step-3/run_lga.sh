cd results
for j in $(ls *pdb);do
		/misc/chaoyi/utils/runlga.mol_mol.pl $j ../../../../target09.resaved.pdb -3 -sia
done
