
options(error=dump.frames)

library(ggplot2)

ld <- list.files(pattern="^crossval[[:digit:]]+$", include.dirs=TRUE)

scale_x_log2 <- function(...)
{
   scale_x_continuous(..., trans=scales::log2_trans())
}

res <- lapply(ld, function(d) {
   cat(d, "\n")
   
   res.lasso <- res.uni <- res.sec <- NULL

   f <- sprintf("%s/crossval.RData", d)
   cat("loading", f, "\n")
   if(file.exists(f)) {
      load(f)
      res.lasso <- cv.d[, 2:3]
      res.lasso$Method <- "Lasso/Elnet"
   }
   
   f <- sprintf("%s/univar.RData", d)
   cat("loading", f, "\n")
   if(file.exists(f)) {
      load(f)
      res.uni <- res
      colnames(res.uni) <- c("Measure", "NonZero")
      res.uni$Method <- "Univariable"
   }
   
   f <- sprintf("%s/secondstage.RData", d)
   cat("loading", f, "\n")
   if(file.exists(f)) {
      load(f)
      res.sec <- res
      colnames(res.sec) <- c("Measure", "NonZero")
      res.sec$Method <- "SecondStage"
   }
   
   res <- rbind(res.lasso, res.uni, res.sec)
   res$Method <- factor(res$Method)
   
   res <- res[res$NonZero > 0, ]
   res
})

res <- do.call(rbind, res)

g <- ggplot(res, aes(x=NonZero, y=Measure, colour=Method))
g <- g + geom_point(alpha=0.5)
g <- g + scale_x_log2()
g <- g + scale_y_continuous(expression(R^2))
g <- g + theme_bw()
g <- g + stat_smooth(method="loess", size=2)
g <- g + coord_cartesian(ylim=c(-0.05, 0.1))

pdf("R2_all.pdf", width=12)
print(g)
dev.off()

