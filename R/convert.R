

simulate <- function(n=5000, p=1e5, beta=rnorm(p + 1),
      xfile="out.bin", yfile="y.csv")
{
   y <- numeric(n)
   
   f <- file(xfile, "wb")
   for(i in 1:n)
   {
      cat(i, "\r")
      x <- rnorm(p)
      y[i] <- ifelse(runif(1) <= drop(plogis(c(1, x) %*% beta)), 1, 0)
      writeBin(y[i], f)
      writeBin(x, f)
   }
   close(f)
   cat("\n")
   
   write.table(y, file=yfile, sep=",", row.names=FALSE, col.names=FALSE)
}

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
ped2bin <- function(pedfile, xfile, yfile, levl=c("A", "G", "C", "T"))
{
   l <- factor(levl)
   m <- contr.treatment(l)

   fin <- file(pedfile, "rt")
   outx <- file(xfile, "wb")
   #outy <- file("
   y <- ""

   while(TRUE)
   {
      r <- readLines(fin, n=1)
      r <- strsplit(r, split=" ")[[1]]
      y <- c(y, r[6])
      seq <- r[-(1:6)]
      table(seq)

      
   }
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
      y <- as.numeric(readLines(yin, n=1))
      if(length(y) == 0)
	 break
      r <- readLines(fin, n=1)
      r <- strsplit(r, split=sep)[[1]]
      x <- c(y, as.numeric(t(cont[r,])))
      cat(i, "y=", y, "x[1:10]:", x[1:10], "\n")
      writeBin(x, con=out)
      i <- i + 1
   }
   cat("\n")

   close(out)
   close(fin)
   close(yin)
}

