library(ggplot2)

dir <- "~/Software/hapgen_1.3"
#roots <- c("sim8", "sim7", "sim6")
roots <- c("sim7")
#resdirs <- c("./", "/mnt", "/mnt")
resdirs <- c("./")
#nums <- c(30000, 10000, 4000)
nums <- 1000
p <- 185805
legend <- sprintf(
   "%s/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt",
   dir)
#nexpers <- c(20, 30, 50)
nexpers <- c(10)
maxfits <- c(20)

runperf <- function(f)
{
   capture.output({s <- system(
	 sprintf("~/Software/perf.src/perf -APR -ROC < %s", f),
	 intern=TRUE)})
   p <- as.numeric(sapply(strsplit(s, "[[:space:]]+"), function(x) x[2]))
   names(p) <- sapply(strsplit(s, "[[:space:]]+"), function(x) x[1])
   p
}


for(k in seq(along=roots))
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
   
   # Analyse CD results
   res.cd <- lapply(seq(along=exper), function(k) {
      ex <- exper[[k]]
      dir <- sprintf("%s/results", ex)
   
      # Find all files and sort by numerical ordering
      files <- list.files(pattern="^beta\\.csv\\.",
   	 path=dir, full.names=TRUE)
      if(length(files) == 0)
         stop("no files found")
      id <- sapply(strsplit(files, "\\."), tail, n=1)
      files <- files[order(as.numeric(id))][1:maxfits]
   
      b.cd <- sapply(files, function(f) {
         read.csv(f, header=FALSE)[-1,1]
      })
      
      df.cd <- apply(b.cd, 2, function(x) sum(x != 0))
   
      mes <- sapply(1:ncol(b.cd), function(i) {
         f <- sprintf("%s/b.cd.%s", dir, i - 1)
         if(!file.exists(f)) 
   	 write.table(cbind(ysnp[[k]], abs(sign(b.cd[,i]))),
   	    col.names=FALSE, row.names=FALSE, sep="\t",
   	    file=f)
         cat(f, "\n")
	 gc()
         runperf(f)
      })
   
      list(measure=mes, df=df.cd, nsim=length(id))
   })
   
   m.cd <- data.frame(
      do.call("rbind", lapply(res.cd, function(r) t(do.call("rbind", r))))
   )
   m.cd$Sim <- factor(
      rep(1:length(res.cd), sapply(res.cd, function(x) length(x$df)))
   )
   
   # Analyse plink results
   res.pl <- lapply(seq(along=exper), function(k) {
      ex <- exper[[k]]
      dir <- sprintf("%s/results", ex)
      d <- read.csv(sprintf("%s/plink.assoc.logistic", dir), sep="")
       
      stats <- cbind(coef=d$STAT, logpval=-log10(d$P))
   
      mes <- sapply(1:ncol(stats), function(i) {
         f <- sprintf("%s/b.cd.plink.%s", dir, colnames(stats)[i])
         if(!file.exists(f))
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
   
   # Both -log10(pval) and STAT yield same AROC/APRC, so take one
   m.pl <- data.frame(t(sapply(res.pl, function(x) x[[1]][,1])))
   # fake DF, to make the point plot nicely
   pl.df <- 300
   m.pl$df <- pl.df
   m.pl$Sim <- 1
   m.pl$nsim <- 1
   m.pl$Method <- "logistic"
   
   m.cd$Method <- "lasso"
   
   # cutoff, don't show the long tail
   m.cd.2 <- m.cd[m.cd$df > 0 & m.cd$df < 256 ,]
   m.comb <- rbind(m.cd.2, m.pl)
   m.comb$Method <- factor(m.comb$Method)
   
   maxdf <- round(max(m.comb$df / 10))
   m.comb$df_bin <- cut(m.comb$df, breaks=(0:maxdf) * 10,
         labels=(1:maxdf - 1) * 10 + 5)
   l <- levels(m.comb$df_bin)
   l[length(l)] <- ""
   levels(m.comb$df_bin) <- l
   
   b <- 2^sort(unique(round(log2(m.comb$df[m.comb$df > 0]))))
   g <- ggplot(m.comb, aes(x=df_bin, y=APR, shape=Method))
   g <- g + geom_boxplot(outlier.size=0)
   g <- g + geom_point(size=3)
   g <- g + ylim(0, 0.42)
   
   pdf(sprintf("%s/results_%s.pdf", resdirs[k], root), width=11)
   print(g)
   dev.off()
   
   save.image(file=sprintf("%s/eval_%s.RData", resdirs[k], root))
}

