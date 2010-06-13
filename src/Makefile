
.PHONY: all

all: sgd pcor

# what about *.h?
#
CCFLAGS = -Wall -g -pg -std=gnu99
LIBRARIES = -lm

sgd: common.c gmatrix.c loss.c evaluation.c util.c sgd.c cd.c gd.c sgd_main.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o sgd

pcor: common.c gmatrix.c loss.c evaluation.c util.c sgd.c cd.c gd.c pcor.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f sgd pcor
# DO NOT DELETE
