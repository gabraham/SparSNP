
PLINK="p-link"

DIR_STEM="sim7."
HAPMAP="$PWD/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt.map"

results="results"
for i in $(seq 31 40);
do
   dir="$DIR_STEM""$i"
   while ! [[ -d "$dir" && -a "$dir/sim.bed" ]];
   do
      echo "Dir $dir doesn't exist yet, sleeping..."
      sleep 660
   done

   dir_res="$dir/$results"
   if ! [ -d "$dir_res" ]; then
      mkdir "$dir_res"
   fi
   pushd "$dir/$results"
   if ! [ -a "plink.assoc.logistic" ];
   then
      eval "$PLINK --bfile ../sim --logistic"
   fi
   popd
done

