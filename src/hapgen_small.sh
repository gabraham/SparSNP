source runhapgen.sh

RR1=1.5
RR2=2.25

for ((J=1 ; J<=1; J++));
do
   DIR="sim2.$J"
   if ! [ -d "$DIR" ]; then
      mkdir $DIR
   fi  
   prefix="sim"
   N=250
   K=2
   HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased.100
   LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt.100
   CLEAN=1
   if ! [ -a "$DIR/sim.bin" ];
   then
      simulate $DIR $prefix $N $K $HAPLO $LEGEND $RR1 $RR2 $CLEAN
   else
      echo "Skipping $DIR"
   fi  
   echo "exit:" $J $?
done

