
.PHONY: all

CC=gcc

# Mac 10.6 differences:
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

CCFLAGS = -std=gnu99 -Wall -ggdb3 \
	  -O2 \
	  -msse \
	  -ftree-vectorize \
	  -fno-omit-frame-pointer \
	  -funroll-loops -Winline \
	  -ffast-math \
	  -fstrict-aliasing \
	  -D_FILE_OFFSET_BITS=$(LONG_BIT) \
	  $(OS)

targets = cd scale transpose hapgen2bin hapgen2ped \
	  cbind split unpack univariable subsample
all: $(targets)

debug: CCFLAGS = -Wall \
   		 -ggdb3 \
		 -std=gnu99 \
		 -D_FILE_OFFSET_BITS=$(LONG_BIT) \
		 $(OS)

debug: $(targets)

static: CCFLAGS += -static

static: $(targets)       

LIBRARIES = -llapack -lblas -lgfortran -lm

cd: common.c coder.c ind.c gmatrix.c loss.c util.c options.c main.c cd.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o cd

scale: common.c coder.c ind.c gmatrix.c cd.c scale.c util.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o scale

transpose: common.c coder.c transpose.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o transpose

hapgen2bin: common.c util.c coder.c hapgen2bin.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o hapgen2bin

hapgen2ped: common.c hapgen2ped.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o hapgen2ped

cbind: common.c cbind.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o cbind

split: common.c util.c ind.c split.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o split

unpack: common.c coder.c ind.c gmatrix.c unpack.c util.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o unpack

univariable: common.c coder.c ind.c gmatrix.c loss.c util.c \
	     options.c cd.c svd.c matrix.c thin.c \
	     multivariable.c univariable.c
	$(CC) $(CCFLAGS) -llapack -lblas -lm $^ $(LIBRARIES) -o univariable \
	$(MULDEF)

subsample: common.c util.c coder.c ind.c gmatrix.c subsample.c
	$(CC) $(CCFLAGS) $^ $(LIBRARIES) -o subsample

clean:
	/bin/rm -f $(targets)

