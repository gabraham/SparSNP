PLINK="p-link"

DIR_STEM="sim6."
HAPMAP="$PWD/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt.map"

results="results"
for i in $(seq 1 10);
do
   dir="$DIR_STEM""$i"
   dir_res="$dir/$results"
   if ! [ -d "$dir_res" ]; then
      mkdir "$dir_res"
   fi
   pushd "$dir/$results"
   $PLINK --ped ../sim.ped --map $HAPMAP --logistic
   popd
done

