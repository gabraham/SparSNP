source runhapgen.sh

RR1=1.5
RR2=2.25

for ((J=31 ; J<=40; J++));
do
   DIR="sim7.$J"
   if ! [ -d "$DIR" ]; then
      mkdir $DIR
   fi  
   prefix="sim"
   N=5000
   K=20
   HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased
   LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt
   CLEAN=1
   if ! [ -a "$DIR/sim.bin" ];
   then
      simulate $DIR $prefix $N $K $HAPLO $LEGEND $RR1 $RR2 $CLEAN
   else
      echo "Skipping $DIR"
   fi  
   echo "exit:" $J $?
done

