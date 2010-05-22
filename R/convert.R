

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
      x <- c(y, as.integer(t(cont[r,])))
      cat(i, "y=", y, "x[1:10]:", x[1:10], "\n")
      writeBin(x, con=out)
      i <- i + 1
   }
   cat("\n")

   cat("Number of variables:", length(x) - 1, "\n")

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

