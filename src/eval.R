library(ggplot2)
library(gdata)

dir <- "~/Software/hapgen_1.3"
roots <- c("sim8", "sim7", "sim6")
resdirs <- c("./", "./", "./")
nums <- c(30000, 10000, 1000)
resdirs <- c("./", "./", "./")
nums <- c(30000, 10000, 1000)
p <- 185805
legend <- sprintf(
   "%s/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt",
   dir)
#nexpers <- c(25, 30, 40)
nexpers <- c(2, 2, 2)
resultsdir.cd <- "results4"
resultsdir.plink <- "results"

runperf <- function(f)
{
   capture.output({s <- system(
	 sprintf("~/Software/perf.src/perf -APR -ROC < %s", f),
	 intern=TRUE)})
   p <- as.numeric(sapply(strsplit(s, "[[:space:]]+"), function(x) x[2]))
   names(p) <- sapply(strsplit(s, "[[:space:]]+"), function(x) x[1])
   p
}

condsign <- function(x, s)
{
   if(s) {
      sign(x)
   } else {
      x
   }
}

analyse <- function(k) 
{
   root <- roots[k]
   n <- nums[k]
   nexper <- nexpers[k]
   exper <- paste(resdirs[k], paste(root, ".", 1:nexper, sep=""), sep="/")
   
   # Get indices of causal SNPs
   cat("Getting true SNP loci ... ")
   snplist <- read.csv(legend, sep="\t")
   ysnp <- lapply(exper, function(ex) {
      cat(ex, "\n")
      snps <- as.character(read.csv(sprintf("%s/loci.txt", ex),
            header=FALSE)[,1])
      pos <- sapply(paste("^", snps, "$", sep=""), grep,
            x=as.character(snplist$position))
      as.numeric((1:p) %in% pos)
   })
   cat("done\n")

    # Analyse plink results
   res.pl <- lapply(seq(along=exper), function(k) {
      ex <- exper[[k]]
      dir <- sprintf("%s/%s", ex, resultsdir.plink)
      f <- sprintf("%s/plink.assoc.logistic", dir)
      d <- read.csv(f, sep="")
      cat(f, "\n")
       
      stats <- cbind(coef=d$STAT, logpval=-log10(d$P))
   
      mes <- sapply(1:ncol(stats), function(i) {
         f <- sprintf("%s/b.cd.plink.%s", dir, colnames(stats)[i])
         #if(!file.exists(f))
   	 write.table(cbind(ysnp[[k]], abs(stats[,i])),
   	    col.names=FALSE, row.names=FALSE, sep="\t",
   	    file=f)
         cat(f, "\n")
         gc()
         runperf(f)
      })

      colnames(mes) <- colnames(stats)
      list(measure=mes)
   })

   # Analyse CD results
   res.cd <- lapply(c(sign=TRUE, nosign=FALSE), function(dosign) {
      lapply(seq(along=exper), function(k) {
         ex <- exper[[k]]
         dir <- sprintf("%s/%s", ex, resultsdir.cd)
      
         # Find all files and sort by numerical ordering
         files <- list.files(pattern="^beta\\.csv\\.",
	       path=dir, full.names=TRUE)
         if(length(files) == 0)
            stop("no files found")
         id <- sapply(strsplit(files, "\\."), tail, n=1)
         files <- files[order(as.numeric(id))]
      
         b.cd <- try(sapply(files, function(f) {
            read.csv(f, header=FALSE)[-1,1]
         }))

         df.cd <- apply(b.cd, 2, function(x) sum(x != 0))
      
         mes <- sapply(1:ncol(b.cd), function(i) {
            f <- sprintf("%s/b.cd.%s", dir, i - 1)
	    write.table(cbind(ysnp[[k]], abs(condsign(b.cd[,i], dosign))),
	       col.names=FALSE, row.names=FALSE, sep="\t", file=f)
            cat(f, "\n")
            gc()
            runperf(f)
         })
      
         list(measure=mes, df=df.cd, nsim=length(id))
      })
   })
   
   m.cd.all <- lapply(res.cd, function(rr) {
      d <- data.frame(
	 do.call("rbind", lapply(rr, function(r) t(do.call("rbind", r)))))
      d$Sim <- factor(
	 rep(1:length(rr), sapply(rr, function(x) length(x$df))))
      d$Method <- "lasso"
      d
   })
   
   # Both -log10(pval) and STAT yield same AROC/APRC, so take one
   m.pl <- data.frame(t(sapply(res.pl, function(x) x[[1]][,1])))

   list(m.cd.all=m.cd.all, m.pl=m.pl)
}

lapply(seq(along=roots), function(k) {
   res <- analyse(k)
   save(res, file=sprintf("%s/eval_%s.RData", resdirs[k], roots[k]))
})


