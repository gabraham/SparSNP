

# Create a design matrix for discrete inputs
design.matrix <- function(xfile.in, xfile.out, n, p)
{
   fin <- file(xfile.in, "rb")
   fout <- file(xfile.out, "wb")

   for(i in 1:n)
   {
      #x <- readBin(fin, what="numeric"
   }
}

# Convert plink PED files to binary format
# Creates a design matrix
#ped2bin <- function(pedfile, xfile, yfile, levl=c("A", "G", "C", "T"))
#{
#   l <- factor(levl)
#   m <- contr.treatment(l)
#
#   fin <- file(pedfile, "rt")
#   outx <- file(xfile, "wb")
#   #outy <- file("
#   y <- ""
#
#   while(TRUE)
#   {
#      r <- readLines(fin, n=1)
#      r <- strsplit(r, split=" ")[[1]]
#      y <- c(y, r[6])
#      seq <- r[-(1:6)]
#      table(seq)
#   }
#}

simulate <- function(n=1000, p=100, beta=rnorm(p + 1),
      noise=rnorm(n), outfile="sim.bin", center=TRUE, scale=TRUE)
{
   x <- scale(matrix(rnorm(n * p), n, p), center=center, scale=scale)
   y <- ifelse(runif(n) <= plogis(cbind(1, x) %*% beta + noise), 1, 0)

   out <- file(outfile, "wb")
   for(i in 1:n)
   {
      writeBin(object=c(y[i], as.numeric(x[i,])), con=out)
   }

   close(out)

   invisible(list(x=x, y=y, beta=beta))
}

# Only makes sense for text files
rowcount <- function(fname, blocksize=512)
{
   f <- file(fname, "rt")

   s <- 0
   # Count number of newlines
   while(length(l <- readLines(f, n=blocksize)) > 0)
      s <- s + length(l)

   close(f)
   s
}

# Packs the binary design matrix into an int
pack <- function(x, k=32)
{
   n <- ceiling(length(x) / k)
   p <- integer(n)
   for(i in 1:n)
      p[i] <- as.integer(sum(x[1:k + (i-1) * k] * 2^(1:k-1), na.rm=TRUE))

   p
}

char2binary <- function(x)
{
   x[x == "0"] <- "\u0000"
   x[x == "1"] <- "\u0001"
   x
}

# Read the ASCII HapGen format and convert to a binary format.
#
# Our format is: y_1 x_11 x_12 ... x_1p
#                y_i x_i1 x_i2 ... x_ip
#                y_n x_n1 x_n2 ... x_np
# 
# The genotype calls 0,1,2 are encoded in a binary design matrix. Each SNP
# uses 2 bits. The data is saved as a character type to save space (only one
# byte per characterr), and the binary SNP is encoded so that its machine
# binary representation (not the same as the binary representation in R which
# is a vector of ints) is correctly read by C code as 0 and 1. Otherwise,
# writing a character "0" would be interpreted as 48 in decimal instead of 0,
# and same for "1" which would be interpreted as 49.
hapgen2bin <- function(hgfile, hgyfile, outfile, alleles=c(0, 1, 2), sep=" ")
{
   fin <- file(hgfile, "rt")
   out <- file(outfile, "wb")
   yin <- file(hgyfile, "rt")

   cont <- contr.treatment(factor(alleles))

   i <- 1
   while(TRUE)
   {
      y <- as.integer(readLines(yin, n=1))
      if(length(y) == 0)
	 break
      r <- readLines(fin, n=1)
      r <- strsplit(r, split=sep)[[1]]
      f <- factor(r, levels=alleles)

      # Equivalent to calling model.matrix on the entire matrix, except for
      # the intercept that we don't include here but is added later in SGD.
      bin <- as.integer(t(cont[r,])) 
      x <- as.character(c(y, bin))
      x <- char2binary(x)
      m <- min(length(x), 10) # only show the first few values
      cat(i, "y=", y, "x:", x[1:m], "...\n")
      writeBin(x, con=out)
      i <- i + 1
   }
   cat("\n")

   #cat("Number of variables:", length(x) - 1, "\n")

   close(out)
   close(fin)
   close(yin)
}

bin2smidas <- function(infile, outfile="smidas", n, p,
      type=c("numeric", "integer"))
{
   type <- match.arg(type)

   fin <- file(infile, "rb")
   fout.x <- sprintf("%s.x", outfile)
   fout.y <- sprintf("%s.y", outfile)

   cat(sprintf("%d %d\n", n, p), file=fout.x)

   for(i in 1:n)
   {
      x <- readBin(fin, what=type, n=p + 1)
      cat(p, paste(1:p - 1, x[-1]), "\n", file=fout.x, append=TRUE)
      cat(x[1], "\n", file=fout.y, append=TRUE)
   }

   close(fin)
}

