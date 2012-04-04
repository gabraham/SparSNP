#!/bin/bash

if [ -z "$1" ];
then
   echo "root name missing"
   exit 1
fi

ROOT=$1

NREPS=10

for((i=1;i<=$NREPS;i++))
do

mkdir rep$i
pushd rep$i

cat > run.qsub <<EOF
#!/bin/bash

#PBS -N $(basename $ROOT)_nofailed-polygenic-$i
#PBS -l nodes=1
#PBS -l pvmem=4G
#PBS -m ae
#PBS -l walltime=1:00:00
# PBS -A VR0126

cd \$PBS_O_WORKDIR

awk 'BEGIN{ srand() }
{
   if(rand() < 0.9) {
      print \$1,\$2
   } else {
      print \$6 > "test.pheno"
   }
}' ../$ROOT.fam > training.txt

plink --noweb --bfile ../$ROOT --keep training.txt \
   --logistic --out training

plink --noweb --bfile ../$ROOT --keep training.txt \
   --model --model-trend \
   --out training

awk '/TEST|TREND/{print \$10}' training.model > trend.txt

for p in 0.8 0.5 0.1 0.05 0.01 0.001 0.0001 0.00001 1e-6 1e-7 1e-8 1e-9
do
   awk -vp=\$p '{
      getline v < "trend.txt"
      if(NR > 1 && v < p) {
	 print \$2,\$4,log(\$10)
      }
   }' training.assoc > score.\$p.txt

   plink --noweb --bfile ../$ROOT \
      --remove training.txt --score score.\$p.txt \
      --out testing.\$p

   Rscript ../eval.R testing.\$p.profile | grep -v WARN > auc.\$p.txt 
done

EOF

popd

done
