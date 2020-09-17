cd results/LDA_RESULTS/

rm -f GDT_ST.dat GDT_HA.dat
for i in {10..129}; do
    cat GDT_ST-${i}.dat >> ../../GDT_ST.dat
done

	
