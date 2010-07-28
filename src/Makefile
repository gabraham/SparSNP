
.PHONY: all


targets = cd scale cd_double tabulate
all: $(targets)

# what about *.h?
#
#CCFLAGS = -Wall -pg -g3 -std=gnu99
CCFLAGS = -Wall -g3 -O2 -ftree-vectorize -fno-omit-frame-pointer -std=gnu99
LIBRARIES = -lm

cd: common.c tabulation.c gmatrix.c loss.c util.c main.c cd.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o cd

cd_double: common.c tabulation.c gmatrix.c loss.c util.c main.c cd.c
	gcc -DDTYPE=double $(CCFLAGS) $(LIBRARIES) $^ -o cd_double

scale: common.c tabulation.c gmatrix.c loss.c util.c cd.c scale.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o scale

tabulate: common.c tabulation.c tabulate.c gmatrix.c loss.c util.c cd.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o tabulate

#pcor: common.c gmatrix.c loss.c util.c cd.c pcor.c
#	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f $(targets)
