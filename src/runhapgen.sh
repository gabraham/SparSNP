################################################################################
#
# Simulate SNP data
#
################################################################################
set -u
set -e
shopt -s extglob

TMPDIR=.

PLINK="p-link"
TRANSPOSE="~/Code/cd/src/transpose"
HAPGEN2BIN="~/Code/cd/src/hapgen2bin"
HAPGEN2PED="~/Code/cd/src/hapgen2ped"

# Cut the HapMap data into $num$ regions in different files
function hapmapcut {
   local legend=$1
   local haplo=$2
   local num=$3
   local cutfile=$4
   local rscript=".Rscript.R"
   local w=`cat $1 | wc -l`

   cat > $rscript <<EOF
   
   n <- $num
   w <- $w - 1
   leg <- "$legend"
   hap <- "$haplo"

   f.leg <- file(leg, "rt")
   f.hap <- file(hap, "rt")

   files.leg <- lapply(1:n, function(i) {
      file(sprintf("%s_%s", leg, i), open="w+")
   })

   files.hap <- lapply(1:n, function(i) {
      file(sprintf("%s_%s", hap, i), open="w+")
   })

   # Approximate split into n blocks
   s <- 1
   while(length(unique(s)) < n)
      s <- sort(sample(n, size=w, replace=TRUE))
   write.table(s, file="$cutfile", col.names=FALSE, row.names=FALSE)

   # First do legend files
   # ignore header
   hd <- readLines(f.leg, n=1)
   for(i in 1:n)
      writeLines(hd, con=files.leg[[i]])

   for(i in 1:w)
   {
      r <- readLines(f.leg, n=1)
      cat(summary(files.leg[[s[i]]])\$description, "\n")
      writeLines(r, con=files.leg[[s[i]]])
   }

   close(f.leg)

   for(f in files.leg)
      close(f)

   # Now do haplo files
   
   i <- 1
   while(TRUE)
   {
      r <- readLines(f.hap, n=1)
      if(length(r) == 0)
	 break
      r2 <- strsplit(r, " ")[[1]]
      for(k in 1:n)
      {
	 r3 <- paste(r2[s == k], collapse=" ")
	 cat(summary(files.hap[[k]])\$description, "\n")
	 writeLines(r3, con=files.hap[[k]])
      }
   }

   for(f in files.hap)
      close(f)

   close(f.hap)

EOF

   Rscript $rscript
}

# http://www.perlmonks.org/?node_id=1910
function randomline {
   RANDLINE=`perl -e 'srand; rand($.) < 1 && ($line = $_) while <>;\
   print $line;' $1`
}

function simulate {
   local DIR=$1
   local prefix=$2
   local N=$3
   local K=$4
   local HAPLO=$5
   local LEGEND=$6
   local rr1=$7
   local rr2=$8
   local clean=$9

   local LOCI=$DIR/loci.txt
   local CUTFILE=$DIR/cut.txt
   
   local P=$(($(head -1 $HAPLO | sed 's/ //g' | wc -c) - 1))
   
   if ! [ -d "$DIR" ]; then
      mkdir $DIR
   fi
   
   echo
   echo "####################################"
   echo "Cutting HapMap data"
   echo "####################################"
   
   # Cut hapmap legend files into blocks
   hapmapcut $LEGEND $HAPLO $K $CUTFILE

   # Convert hapmap legend file to plink MAP file
   if ! [ -a "$LEGEND.map" ];
   then
      local RSCRIPT=".conv.R"
      cat > $RSCRIPT <<EOF
   source("~/Code/cd/R/convert.R")
   hapmap2map("$LEGEND", "$LEGEND.map")
EOF
      Rscript $RSCRIPT
   fi


   echo
   echo "####################################"
   echo "Simulating genotypes"
   echo "####################################"
   
   if [ $clean == 1 ];
   then
      /bin/rm -f $LOCI
   fi
   
   # Simulate genotypes using each of the causal SNPs
   for ((i = 1 ; i <= $K ; i++));
   do
      local f="$LEGEND""_$i"
      local w=$(cat $f | wc -l)
      tail -n $((w-1)) $f > "$f"".tmp"
      randomline "$f"".tmp"
      SNP=`echo $RANDLINE | cut -f 2 -d ' '`
      echo $SNP >> $LOCI
   
      CMD="./hapgen -h $HAPLO""_$i -l $LEGEND""_$i \
      -o "$DIR/sim$i" -n $N $N -gen -rr $rr1 $rr2  -dl $SNP"
      echo $CMD
      set +e # hapgen returns 1 on exit
      eval $CMD
      set -e
   
      /bin/rm -f "$f"".tmp"
   done
   
   # Concatenate the hapgen files *column-wise* (same as cbind in R)
   #
   # This depends on hapgen always generating the same response classes
   # for the same samples (which it does)
   #
   RSCRIPT=".Rscript.R"
   cat > $RSCRIPT <<EOF
   sp <- $K
   
   fout <- file("$DIR/sim.all.g", "wt")
   nrow <- 2 * $N

   files <- lapply(1:sp, function(i) {
      file(sprintf("$DIR/sim%s.all.g", i), open="r")
   })

   for(i in 1:nrow)
   {
      cat("row:", i, "\n")
      r <- lapply(files, readLines, n=1)
      r2 <- paste(r, collapse=" ")
      r3 <- gsub("  ", " ", r2)
      writeLines(r3, con=fout)
   }

   close(fout)
EOF

   Rscript $RSCRIPT
   
   # See previous comment re column wise
   /bin/cp $DIR/sim1.y $DIR/sim.y
   if [ $clean == 1 ];
   then
      /bin/rm -rf $DIR/sim+([0-9]).all.g
   fi
   
   echo "####################################"
   echo "Postprocessing genotypes"
   echo "####################################"
   
   # For coordinate descent 
   #convert $DIR $prefix $((N*2))
   cmd="$HAPGEN2BIN -finx $DIR/sim.all.g -finy $DIR/sim.y \
   -fout $DIR/sim.bin -n $((N*2)) -p $P"
   echo $cmd
   eval $cmd

   #echo "####################################"
   #echo "Transposing ..."
   #eval "$TRANSPOSE" -fin "$DIR/sim.bin" -fout "$DIR/sim.bin.t" \
   #-n $((N*2)) -p $((P+1)) 
   #if [ $clean == 1 ];
   #then
   #   /bin/rm "$DIR/sim.bin"
   #fi
   #echo "####################################"

   echo "####################################"
   echo "Converting to plink PED format"
   echo "####################################"

   # For plink, text ped format
#   cat > $RSCRIPT <<EOF
#   source("~/Code/cd/R/convert.R")
#   hapgen2ped("$DIR/sim.all.g", "$DIR/sim.y", "$DIR/sim.ped")
#EOF
#   Rscript $RSCRIPT
   eval "$HAPGEN2PED" -finx "$DIR/sim.all.g" -finy "$DIR/sim.y" \
   -fout "$DIR/sim.ped" -n $((N*2)) -p $P

   if [ $clean == 1 ];
   then
      /bin/rm $DIR/sim.all.g
   fi

   # plink, binary bed format
   $PLINK --ped "$DIR/sim.ped" --map "$LEGEND.map" --make-bed --out "$DIR/sim"
   if [ $clean == 1 ];
   then
      /bin/rm "./$DIR/sim.ped"
   fi

   echo
   echo "####################################"
   echo "DONE"
   echo "####################################"
}

