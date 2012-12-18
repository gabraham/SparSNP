#!/bin/bash

# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
# All rights reserved.
# 

export ROOT=$1

if ! [ -d discovery ];
then
   mkdir discovery
fi

[[ -z "$DOSPARSNP" ]] && DOSPARSNP=1

for((i=1;i<=30;i++))
do
   export REP_START=$i
   export REP_END=$i

cat > run$i.qsub <<EOF
#!/bin/bash

#PBS -N $(basename $ROOT)-crossval$i
#PBS -l nodes=1
#PBS -l pvmem=8G
#PBS -m ae
#PBS -l walltime=3:00:00:00

cd \$PBS_O_WORKDIR

if [ "$DOSPARSNP" ];
then
   crossval.sh \$ROOT sqrhinge 2>&1 | tee sparsnp.log
fi

if [ "$DOUNIVAR" ];
then
   crossvaluni.sh \$ROOT 2>&1 | tee univar.log
fi

EOF

   if ! [ -d "discovery/crossval$i" ];
   then
      qsub -V run$i.qsub
   fi
   

done

