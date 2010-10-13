source runhapgen.sh

RR1=1.5
RR2=2.25

for ((J=4 ; J<=4; J++));
do
   DIR="sim3.$J"
   if   [ -d "$DIR" ]; then
      mkdir $DIR
   fi  
   prefix="sim"
   N=1000
   K=20
   HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased.1000
   LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt.1000
   CLEAN=0
   if ! [ -a "$DIR/sim.bin" ];
   then
      simulate $DIR $prefix $N $K $HAPLO $LEGEND $RR1 $RR2 $CLEAN
   else
      echo "Skipping $DIR"
   fi  
   echo "exit:" $J $?
done

