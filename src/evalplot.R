
library(ggplot2)
library(gdata)

plot.apr2 <- function(m.cd, m.pl, suf, nz.breaks, maxdf=256)
{
   # cutoff, don't show the long tail
   m.pl$df <- 1
   m.cd.2 <- m.cd[m.cd$df > 0 & m.cd$df < maxdf ,]
   m.comb <- rbind(m.cd.2, m.pl)
   m.comb$Method <- factor(m.comb$Method)
   
   #rmaxdf <- round(max(m.comb$df / 20))
   #m.comb$df_bin <- cut(m.comb$df, breaks=(0:rmaxdf) * 20)
   m.comb$df_bin <- cut(m.comb$df, breaks=nz.breaks, include.lowest=TRUE)
   l <- levels(m.comb$df_bin)
   l <- gsub(",", "--", gsub("\\[|\\(|\\]", "", l))
   d <- diff(c(0, nz.breaks))[-1]
   w <- which(d == 1)
   l[w] <- nz.breaks[w]
   levels(m.comb$df_bin) <- l
   m.comb$df_bin <- drop.levels(m.comb$df_bin, reorder=FALSE)
   levels(m.comb$Method) <- c("lasso", "logis")

   g <- ggplot(m.comb, aes(x=df_bin, y=APR, shape=Method))
   g <- g + geom_boxplot(outlier.shape=as.numeric(NA),
	 colour="darkgray", size=1.5)
   g <- g + geom_point(size=3, position=position_jitter(width=0.07))
   g <- g + ylim(0, 0.42)
   g <- g + scale_shape_manual(values=c(1, 2))
   g <- g + xlab("# SNPs with non-zero effects") + ylab("APRC")
   g <- g + opts(
	 axis.text.x=theme_text(size=20, angle=90, hjust=1),
	 axis.text.y=theme_text(size=23, angle=0, hjust=0),
	 axis.title.x=theme_text(size=23),
	 axis.title.y=theme_text(size=23, angle=90),
	 axis.ticks=NULL,
	 strip.text.x=theme_text(size=20),
	 #strip.text.y=theme_text(size=23),
	 legend.position="none")
   g <- g + facet_grid(~ Method, scales="free", space="free")
   
   pdf(sprintf("hapgen_APR_%s.pdf", suf), width=11)
   print(g)
   dev.off()
}


legend <- "~/Software/hapgen_1.3/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt"
nums <- c(30000, 10000, 1000)
resultsdir.cd <- "results4"
roots <- c("sim8", "sim7", "sim6")
nz.breaks <- c(1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 512)

lapply(seq(along=roots), function(i) {
   load(sprintf("eval_%s.RData", roots[i]))
   rm(maxdf, suf)

   # Use the second version, nosign
   plot.apr2(m.cd.all[[2]], m.pl, roots[i], nz.breaks)
})

