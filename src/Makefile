
.PHONY: all

all: sgd

# what about *.h?
#
CCFLAGS = -Wall -g -pg -std=gnu99
LIBRARIES = -lm

sgd: common.c gmatrix.c scale.c loss.c evaluation.c util.c sgd.c sgd_main.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o sgd

pcor: common.c gmatrix.c scale.c loss.c evaluation.c util.c sgd.c pcor.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f sgd #sgd_cont sgd_disc pcor
# DO NOT DELETE
