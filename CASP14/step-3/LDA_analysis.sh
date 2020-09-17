cd results/LDA_RESULTS
for i in {10..129}; do
	        cat target09-${i}.pdb.*.pdb.res | grep "GDT PERCENT_AT" | awk '{ V=($4+$6+$10+$18)/4.0; printf "%6.2f\n",V; }' >> GDT_ST-${i}.dat
done

