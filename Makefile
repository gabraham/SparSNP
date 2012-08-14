
.PHONY: all

#CC=gcc

# Mac OSX differences:
#
# * linker won't accept -allow-multiple-definition, required on Linux 
#   as a hack for xerbla clash between lapack and blas
# * Accelerate.h headers
MULDEF = -Wl,--allow-multiple-definition 
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   OS = -DMACOSX
   MULDEF =
endif

# Allow large-file seeking with fseeko
LONG_BIT = $(shell getconf LONG_BIT)

CFLAGS = -std=gnu99 -Wall -ggdb3 -g3 \
	  -O3 \
	  -msse \
	  -ftree-vectorize \
	  -fno-omit-frame-pointer \
	  -funroll-loops -Winline \
	  -ffast-math \
	  -fstrict-aliasing \
	  -D_FILE_OFFSET_BITS=$(LONG_BIT) \
	  $(OS)

targets = sparsnp scale transpose \
	  cbind makefolds unpack univariable subsample \
	  realpath gennetwork_test

all: $(targets)

debug: CFLAGS = -Wall \
	        -DDEBUG \
		-ggdb3 \
		-std=gnu99 \
		-D_FILE_OFFSET_BITS=$(LONG_BIT) \
		$(OS)

debug: $(targets)


# static linking is NOT supported on OSX
static: CFLAGS += -static

static: $(targets)       

LIBRARIES = -lpthread -llapack -lblas -lm #-lgfortran -lm

sparsnp: common.c coder.c ind.c gmatrix.c link.c util.c options.c \
	 main.c sparsnp.c matrix.c gennetwork.c
	$(CC) $(CFLAGS) $^ $(LIBRARIES) -o sparsnp

scale: common.c coder.c ind.c gmatrix.c sparsnp.c scale.c util.c \
	 matrix.c gennetwork.c
	$(CC) $(CFLAGS) $^ $(LIBRARIES) -o scale

transpose: common.c coder.c transpose.c
	$(CC) $(CFLAGS) $^ $(LIBRARIES) -o transpose

cbind: common.c cbind.c
	$(CC) $(CFLAGS) $^ $(LIBRARIES) -o cbind

makefolds: common.c util.c ind.c makefolds.c matrix.c
	$(CC) $(CFLAGS) $^ $(LIBRARIES) -o makefolds

unpack: common.c coder.c ind.c gmatrix.c unpack.c util.c \
	 matrix.c gennetwork.c
	$(CC) $(CFLAGS) $^ $(LIBRARIES) -o unpack

univariable: common.c coder.c ind.c gmatrix.c link.c util.c \
	     options.c sparsnp.c svd.c matrix.c thin.c \
	     multivariable.c univariable.c gennetwork.c
	$(CC) $(CFLAGS) -llapack -lblas -lm $^ $(LIBRARIES) -o univariable \
	$(MULDEF)

subsample: common.c util.c coder.c ind.c gmatrix.c subsample.c \
	 matrix.c gennetwork.c
	$(CC) $(CFLAGS) $^ $(LIBRARIES) -o subsample

realpath: realpath.c
	$(CC) $(CFLAGS) $^ -o realpath

gennetwork_test: gennetwork.c gennetwork_test.c matrix.c util.c
	$(CC) $(CFLAGS) -ggdb3 $^ $(LIBRARIES) -o gennetwork_test

clean:
	/bin/rm -f $(targets)

