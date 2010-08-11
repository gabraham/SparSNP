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

# Row major ordering
#Line 1: m d
#Line 2: i1 v1 i2 v2 ...
#Line d+1: i1 v1 i2 v2 ...
function formatsmidas {
   local DIR=$1
   local binfile=$2
   local N=$3
   local P=$4
   local RSCRIPT=formatsmidas.R
   local FILEX=smidas_x.dat
   local FILEY=smidas_y.dat
   echo "formatsmidas: $DIR $binfile $N $P"

   cat > $DIR/$RSCRIPT <<EOF
   n2 <- $N * 2
   p <- $P
   xin <- file("$DIR/$binfile", "rb")
   cat(n2, p + 1, "\n", file="$DIR/$FILEX")
   cat("", file="$DIR/$FILEY")

   for(i in 1:n2)
   {
      r <- as.numeric(readBin(con=xin, what="raw", n=p + 1))
      cat(p + 1, paste(0:p, c(1, r[-1])), "\n",
	    file="$DIR/$FILEX", append=TRUE)
      cat(ifelse(r[1] == 1, 1, -1), "\n",
	    file="$DIR/$FILEY", append=TRUE)
   }
   close(xin)

EOF

   Rscript $DIR/$RSCRIPT
   /bin/rm $DIR/$RSCRIPT

}

# Column-major ordering
function formatscd {
   local DIR=$1
   local binfile=$2
   local N=$3
   local P=$4
   local RSCRIPT=formatscd.R
   local FILEX=scd_x.dat
   local FILEY=scd_y.dat
   echo "formatscd: $DIR $binfile $N $P"

   cat > $DIR/$RSCRIPT <<EOF
   n2 <- $N * 2
   p <- $P
   xin <- file("$DIR/$binfile", "rb")
   cat(n2, p + 1, "\n", file="$DIR/$FILEX")
   cat("", file="$DIR/$FILEY")

   cat(paste(1:n2 - 1, rep(1, n2)), "\n", file="$DIR/$FILEX", append=TRUE)

   for(j in 1:p)
   {
      for(i in 1:n2)
      {
         r <- as.numeric(readBin(con=xin, what="raw", n=p + 1))

         cat(i - 1, r[j + 1], file="$DIR/$FILEX", append=TRUE)
	 if(i < n2) {
	    cat(" ", file="$DIR/$FILEX", append=TRUE)
	 } else {
	    cat("\n", file="$DIR/$FILEX", append=TRUE)
	 }

	 if(j == 1)
	 {
	    cat(r[1], "\n")
	    cat(ifelse(r[1] == 1, 1, -1), "\n",
               file="$DIR/$FILEY", append=TRUE)
	 }
      }
      close(xin)
      xin <- file("$DIR/$binfile", "rb")
   }
   close(xin)

EOF

   Rscript $DIR/$RSCRIPT
   /bin/rm $DIR/$RSCRIPT

}

function formatsvmlight {
   local DIR=$1
   local binfile=$2
   local N=$3
   local P=$4
   local RSCRIPT=formatsvmlight.R
   local FILEX=svmlight_x.dat
   echo "formatsvmlight: $DIR $binfile $N $P"

   cat > $DIR/$RSCRIPT <<EOF
   n2 <- $N * 2
   p <- $P
   xin <- file("$DIR/$binfile", "rb")
   cat("", file="$DIR/$FILEX")

   k <- 1:(p+1)
   for(i in 1:n2)
   {
      r <- as.numeric(readBin(con=xin, what="raw", n=p + 1))
      y <- ifelse(r[1] == 1, 1, -1) 
      v <- c(1, r[-1])
      w <- which(v != 0)

      cat(y, paste(k[w], v[w], sep=":"), "\n",
	    file="$DIR/$FILEX", append=TRUE)
   }
   close(xin)

EOF

   Rscript $DIR/$RSCRIPT
   /bin/rm $DIR/$RSCRIPT

}

function testscd {
   local DIR=testscd
   local RSCRIPT=testscd.R
   local binfile=test.bin
   local n=6
   local p=$((2*6))
   
   if ! [ -d "$DIR" ]; then
      mkdir $DIR
   fi

   cat > $DIR/$RSCRIPT <<EOF
   n <- $n * 2
   p <- $p
   y <- as.raw(rep(c(0, 1), each=n/2))
   x <- t(matrix(as.raw(rep(c(0, 1, 0), each=n, times=4)), n, p))
   f <- file("$DIR/$binfile", "wb")
   for(i in 1:n)
   {
      writeBin(c(y[i], x[i,]), con=f)
   }
   close(f)
EOF
   
   Rscript $DIR/$RSCRIPT

   formatscd $DIR "$binfile" $n $p
}

function convert {
   local DIR=$1
   local prefix=$2
   local N=$3
   local binfile="$prefix.bin"
   local xfile="$prefix.all.g"
   local yfile="$prefix.y"
   local rscript=".convert.R"
   
   cat > $rscript <<EOF
   hgfile <- "$DIR/$xfile"
   hgyfile <- "$DIR/$yfile"
   outfile <- "$DIR/$binfile"
   sep <- " "

   fin <- file(hgfile, "rt")
   yin <- file(hgyfile, "rt")
   out <- file(outfile, "wb")

   i <- 1
   while(TRUE)
   {
      y <- as.integer(readLines(yin, n=1))
      if(length(y) == 0)
	 break
      r <- readLines(fin, n=1)
      r <- strsplit(r, split=sep)[[1]]

      cat("read", length(r), "fields\n")

      x <- as.raw(c(y, r))
      cat(i, "y=", y, length(x), "following: x:", x[1:21], "...\n")
      writeBin(x, con=out)
      flush(out)
      i <- i + 1
   }

   close(fin)
   close(yin)
   warnings()
   
EOF

   Rscript $rscript

   #/bin/rm -f $DIR/$xfileshuf $rscript
 
}

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

function transpose {
   local infile="$1"
   local outfile="$infile.t"
   local n=$2
   local p=$3
   local rscript=.transpose.R

   cat > $rscript <<EOF
   x <- matrix(readBin("$infile", what="raw",
	 n=$n * ($p + 1)), nrow=$n, byrow=TRUE)
   gc()
   writeBin(as.raw(x), con="$outfile")
EOF
   Rscript $rscript
   /bin/rm $rscript
}

function simulate {
   local DIR=$1
   local prefix=$2
   local N=$3
   local K=$4
   local HAPLO=$5
   local LEGEND=$6

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
   
   /bin/rm -f $LOCI
   
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
      -o "$DIR/sim$i" -n $N $N -gen -rr 1.5 2.25  -dl $SNP"
      echo $CMD
      set +e # hapgen returns 1 on exit
      eval $CMD
      set -e
   
      /bin/rm -f "$f"".tmp"
   done
   
   # Concatenate the hapgen files *column-wise*
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
   #for(f in files)
   #   close(f)
EOF

   Rscript $RSCRIPT
   
   # See previous comment
   /bin/cp $DIR/sim1.y $DIR/sim.y
   /bin/rm -rf $DIR/sim+([0-9]).all.g
   
   
   echo "####################################"
   echo "Postprocessing genotypes"
   echo "####################################"
   
   # For coordinate descent 
   convert $DIR $prefix $((N*2))
   transpose "$DIR/sim.bin" $((N*2)) $P 
   /bin/rm "$DIR/sim.bin"

   echo "####################################"
   echo "Converting to plink PED format"
   echo "####################################"

   # For plink, text ped format
   cat > $RSCRIPT <<EOF
   source("~/Code/cd/R/convert.R")
   hapgen2ped("$DIR/sim.all.g", "$DIR/sim.y", "$DIR/sim.ped")
EOF
   Rscript $RSCRIPT

   rm $DIR/sim.all.g

   # plink, binary bed format
   $PLINK --ped "$DIR/sim.ped" --map "$LEGEND.map" --make-bed --out "$DIR/sim"
   /bin/rm "./$DIR/sim.ped"

   echo
   echo "####################################"
   echo "DONE"
   echo "####################################"
}

#################################################################################


# HapMap data, one strong SNP

#DIR=sim1
#prefix="sim"
#N=1000 # No. samples in each group
#HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased
#LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt
#SNP=72434
#
#if ! [ -d "$DIR" ]; then
#   mkdir $DIR
#fi
#
#set +e
#./hapgen -h $HAPLO -l $LEGEND \
#-o $DIR/$prefix -n $N $N -gen -rr 1.5 2.25 -dl $SNP
#set -e
#
##convert $DIR $prefix $((2*N))
#
#exit 1


################################################################################

for ((J=1 ; J<=10; J++));
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
   simulate $DIR $prefix $N $K $HAPLO $LEGEND
   echo "exit:" $J $?
done

exit 1

## HapMap data, several strong SNPs
#DIR=sim7
#prefix="sim"
#N=1000
#K=20
#HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased
#LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt
#simulate $DIR $prefix $N $K $HAPLO $LEGEND

##
## Several strong SNPs, lots of weak SNPs
##
#



