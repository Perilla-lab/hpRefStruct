x=46
for i in {10..129} ; do
cat results/score/target$x-$i.sc | grep SCORE | awk '{print $2}' | tail -1 >> energy
done
