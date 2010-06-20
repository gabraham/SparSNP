
.PHONY: all

all: sgd sgd_double #pcor

# what about *.h?
#
CCFLAGS = -Wall -g -pg -std=gnu99
LIBRARIES = -lm

sgd: common.c gmatrix.c loss.c evaluation.c util.c sgd.c cd.c gd.c sgd_main.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o sgd

sgd_double: common.c gmatrix.c loss.c evaluation.c util.c sgd.c cd.c gd.c sgd_main.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o sgd_double

pcor: common.c gmatrix.c loss.c evaluation.c util.c sgd.c cd.c gd.c pcor.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f sgd pcor
# DO NOT DELETE
