x=46

mkdir -p results/score


for i in {10..129}; do
/home_orig/chaoyi/builds/rosetta_src_2020.08.61146_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/score_jd2.default.linuxgccrelease -in:file:s results/target${x}-${i}.pdb -out:file:scorefile results/score/target${x}-${i}.sc -ignore_unrecognized_res
done

