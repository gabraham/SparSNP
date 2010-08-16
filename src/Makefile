
.PHONY: all

CCFLAGS = -std=gnu99 -Wall -ggdb3 -O2 -ftree-vectorize -fno-omit-frame-pointer \
	  -funroll-loops -Winline \
	  -ffast-math -fstrict-aliasing

targets = cd scale transpose
all: $(targets)

debug: CCFLAGS = -Wall -ggdb3 -std=gnu99
debug: $(targets)

LIBRARIES = -lm

cd: common.c tabulation.c gmatrix.c loss.c util.c options.c main.c cd.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o cd

scale: common.c gmatrix.c util.c cd.c scale.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o scale

transpose: common.c transpose.c
	$(CC) $(CCFLAGS) $(LIBRARIES) $^ -o transpose

#pcor: common.c gmatrix.c loss.c util.c cd.c pcor.c
#	$(CC) -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f $(targets)
