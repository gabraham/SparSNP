################################################################################
#
# Simulate SNP data
#
################################################################################
set -u
set -e

TMPDIR=.


# Row major ordering
#Line 1: m d
#Line 2: i1 v1 i2 v2 ...
#Line d+1: i1 v1 i2 v2 ...
function formatsmidas {
   DIR=$1
   binfile=$2
   N=$3
   P=$4
   RSCRIPT=formatsmidas.R
   FILEX=smidas_x.dat
   FILEY=smidas_y.dat
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
   DIR=$1
   binfile=$2
   N=$3
   P=$4
   RSCRIPT=formatscd.R
   FILEX=scd_x.dat
   FILEY=scd_y.dat
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
   DIR=$1
   binfile=$2
   N=$3
   P=$4
   RSCRIPT=formatsvmlight.R
   FILEX=svmlight_x.dat
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
   DIR=testscd
   RSCRIPT=testscd.R
   binfile=test.bin
   n=6
   p=$((2*6))
   
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

function shuffle {
   DIR=$1
   prefix=$2
   N=$3
   binfile="$prefix.bin"
   xfile="$prefix.all.g"
   xfileshuf="$xfile.shuffled"
   yfile="$prefix.y"
   yfileshuf="$yfile.shuffled"
   tmp1=".tmp1"
   tmp2=".tmp2"
   rscript=".shuffle.R"
   
   echo -n "Shuffling ... "
   # Shuffle samples 
   cat > $rscript <<EOF
   n2 <- $N
   s <- sample(n2)
   y <- read.csv("$DIR/$yfile", header=FALSE)[,1]
   write.table(y[order(s)], "$DIR/$yfileshuf", col.names=FALSE, row.names=FALSE)
   f <- file("$DIR/$xfile", "rt")
   tmp <- file("$DIR/$tmp1", "wt")
   for(i in 1:n2)
   {
      r <- readLines(f, n=1)
      writeLines(paste(s[i], r, sep=" "), tmp)
   }
   close(tmp)
   close(f)
EOF
   
   Rscript $rscript
   
   sort -T $TMPDIR -n -k 1,1 $DIR/$tmp1 > $DIR/$tmp2
   /bin/rm $DIR/$tmp1

   cut -f2- -d ' ' $DIR/$tmp2 > $DIR/$xfileshuf
   /bin/rm $DIR/$tmp2 $rscript
   
}

function convert {
   DIR=$1
   prefix=$2
   N=$3
   binfile="$prefix.bin"
   xfile="$prefix.all.g"
   xfileshuf="$xfile.shuffled"
   yfile="$prefix.y"
   yfileshuf="$yfile.shuffled"
   rscript=".convert.R"
   
   cat > $rscript <<EOF
   hgfile <- "$DIR/$xfileshuf"
   hgyfile <- "$DIR/$yfileshuf"
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
   legend=$1
   haplo=$2
   num=$3
   cutfile=$4
   rscript=".Rscript.R"

   w=`cat $1 | wc -l`
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

#################################################################################
## Small example with one strong SNP
#
#DIR="ex"
#prefix="sim"
#
## No. samples in each group
#N=5000
#
#if ! [ -d "$DIR" ]; then
#   mkdir $DIR
#fi
#
#./hapgen -h example/ex.haps \
#-l example/ex.leg \
#-r example/ex.map \
#-o $DIR/$prefix -n $N $N -gen -rr 2 4 -dl 14431347
#
#shuffle $DIR $prefix $((2*N))
#convert $DIR $prefix $((2*N))
#
#exit 1

################################################################################

#testscd

#exit 1

#################################################################################
## Truncated HapMap data, one strong SNP
#
#DIR=sim1
#prefix="sim"
#p=100
#N=1000 # No. samples in each group
#HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased.100
#LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt.100
#SNP=72434
#
#if ! [ -d "$DIR" ]; then
#   mkdir $DIR
#fi
#
#set +e
#./hapgen -h $HAPLO -l $LEGEND \
#-o $DIR/$prefix -n $N $N -gen -rr 4 8 -dl $SNP
#set -e
#
#shuffle $DIR $prefix $((2*N))
#convert $DIR $prefix $((2*N))
#
#exit 1

##formatsmidas $DIR "$prefix.bin" $N $p
##formatscd $DIR "$prefix.bin" $N $p
#
#
################################################################################
# Truncated HapMap data, one strong SNP, different SNP
#
#DIR=sim2
#SNP=555296
#
#if ! [ -d "$DIR" ]; then
#   mkdir $DIR
#fi
#
#./hapgen -h $HAPLO -l $LEGEND \
#-o $DIR/$prefix -n $N $N -gen -rr 4 8 -dl $SNP
#
#shuffleconv $DIR $prefix $N
#
###formatsmidas $DIR "$prefix.bin" $N $p
###formatscd $DIR "$prefix.bin" $N $p
##
#
#exit 1
#
#################################################################################
## HapMap data, one strong SNP
#
#DIR=sim3
#prefix="sim"
#
#if ! [ -d "$DIR" ]; then
#   mkdir $DIR
#fi
#
## No. samples in each group
#N=2500
#p=185805
#HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased
#LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt
#
##./hapgen -h $HAPLO -l $LEGEND \
##-o $DIR/$prefix -n $N $N -gen -rr 4.0 8.0 -dl 555296
#
##shuffleconv $DIR $prefix $N
#
#formatsvmlight $DIR "$prefix.bin" $N $p
#
#exit 1

#################################################################################
# HapMap data, several strong SNPs

DIR=sim5
prefix="sim"
N=500
HAPLO=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd.phased.50000
LEGEND=HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt.50000
LOCI=$DIR/loci.txt
CUTFILE=$DIR/cut.txt

if ! [ -d "$DIR" ]; then
   mkdir $DIR
fi

# Number of causal SNPs
K=10

echo
echo "####################################"
echo "Cutting HapMap data"
echo "####################################"

hapmapcut $LEGEND $HAPLO $K $CUTFILE

<<<<<<< HEAD:src/runhapgen.sh
#exit 1
=======
echo
echo "####################################"
echo "Simulating genotypes"
echo "####################################"
echo

/bin/rm -f $LOCI

# Simulate genotypes using each of the causal SNPs
for ((i = 1 ; i <= $K ; i++));
do
   f="$LEGEND""_$i"
   w=$(cat $f | wc -l)
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


echo "####################################"
echo "Postprocessing genotypes"
echo "####################################"

shuffle $DIR $prefix $((N*2))
convert $DIR $prefix $((N*2))
