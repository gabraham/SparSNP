
.PHONY: all

all: cd cd_double #pcor

# what about *.h?
#
CCFLAGS = -Wall -g3 -pg -std=gnu99
LIBRARIES = -lm

cd: common.c gmatrix.c loss.c evaluation.c util.c cd.c gd.c main.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o cd

cd_double: common.c gmatrix.c loss.c evaluation.c util.c cd.c gd.c main.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o cd_double

pcor: common.c gmatrix.c loss.c evaluation.c util.c cd.c gd.c pcor.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f cd cd_double pcor
# DO NOT DELETE
