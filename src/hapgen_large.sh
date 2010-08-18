source runhapgen.sh

RR1=1.5
RR2=2.25

for ((J=1 ; J<=1; J++));
do
   DIR="sim8.$J"
   if ! [ -d "$DIR" ]; then
      mkdir $DIR
   fi  
   prefix="sim"
   N=25000
   K=20
   HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased
   LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt
   if ! [ -a "$DIR/sim.bin.t" ];
   then
      simulate $DIR $prefix $N $K $HAPLO $LEGEND $RR1 $RR2
   else
      echo "Skipping $DIR"
   fi  
   echo "exit:" $J $?
done
