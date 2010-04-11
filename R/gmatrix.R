setClass("gmatrix",
   representation(env="environment", companions="list"),
   prototype(env=new.env(), companions=list())) 

setClass("gmatrixMem",
   representation(x="matrix", nrow="integer", ncol="integer"),
   prototype(x=matrix(), nrow=0L, ncol=0L),
   contains="gmatrix")

setOldClass("file")

setClass("gmatrixDisk",
   representation(x="character", file="file", ncol="integer"),
   prototype(x=character(), file=file(), ncol=0L),
   contains="gmatrix")

# companions: a list of variabels to be indexed together with x
gmatrixMem <- function(x, nrow=nrow(x), ncol=ncol(x), companions=list())
{
   g <- new("gmatrixMem", x=x, nrow=nrow, ncol=ncol, companions=companions)
   g@env[["i"]] <- 1
   g
}

gmatrixDisk <- function(x, nrow, ncol, companions=list())
{
   g <- new("gmatrixDisk", x=x, nrow=nrow, ncol=ncol, file=file(x, "rb"),
	 companions=companions)
   g@env[["i"]] <- 1
   g
}

setGeneric("nextRow", function(object, ...) standardGeneric("nextRow"))

setMethod("nextRow", signature("gmatrixDisk"),
   function(object, loop=TRUE, companions=TRUE, blocksize=1) {
      dat <- readBin(object@file, what="numeric", n=object@ncol)
      if(length(dat) == 0) {
	 if(loop) {
	    object@file <- file(object@x, "rb")
	    dat <- readBin(object@file, what="numeric", n=object@ncol)
	 } else {
	    return(numeric(0))
	 }
      }
      object@env[["i"]] <- object@env[["i"]] + 1L
      m <- matrix(dat, nrow=min(blocksize, length(dat) / object@ncol),
	 ncol=object@ncol, byrow=TRUE)
      if(companions && length(object@companions) > 0) {
	 list(x=m,
	    companions=lapply(object@companions, function(x) {
	       x[object@env[["i"]] - 1L, , drop=FALSE]
	    }))
      } else m
   }
)

# companions: logical, return companions list or not
setMethod("nextRow", signature("gmatrixMem"),
   function(object, loop=TRUE, companions=TRUE) {
      if(object@env[["i"]] > nrow(object@x)) {
	 if(loop) {
	    object@env[["i"]] <- 1L
	 } else {
	    return(numeric(0))
	 }
      }
      r <- object@x[object@env[["i"]], , drop=FALSE] 
      object@env[["i"]] <- object@env[["i"]] + 1L
      if(companions && length(object@companions) > 0) {
	 list(x=r,
	    companions=lapply(object@companions, function(x) {
	       x[object@env[["i"]] - 1L, , drop=FALSE]
	    }))
      } else r
   }
)

