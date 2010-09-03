
.PHONY: all

#CC=gcc

# Allow large-file seeking with fseeko
LONG_BIT = $(shell getconf LONG_BIT)

CCFLAGS = -std=gnu99 -Wall -ggdb3 \
	  -O2 \
	  -ftree-vectorize \
	  -fno-omit-frame-pointer \
	  -funroll-loops -Winline \
	  -ffast-math \
	  -fstrict-aliasing \
	  -D_FILE_OFFSET_BITS=$(LONG_BIT)

targets = cd scale transpose hapgen2bin hapgen2ped cbind
all: $(targets)

debug: CCFLAGS = -Wall -ggdb3 -std=gnu99 -D_FILE_OFFSET_BITS=$(LONG_BIT)
debug: $(targets)

LIBRARIES = -lm

cd: common.c coder.c gmatrix.c loss.c util.c options.c main.c cd.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o cd

scale: common.c coder.c gmatrix.c util.c cd.c scale.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o scale

transpose: common.c coder.c transpose.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o transpose

hapgen2bin: common.c coder.c hapgen2bin.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o hapgen2bin

hapgen2ped: common.c hapgen2ped.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o hapgen2ped

cbind: common.c cbind.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o cbind

#pcor: common.c gmatrix.c loss.c util.c cd.c pcor.c
#	$(CC) -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f $(targets)
