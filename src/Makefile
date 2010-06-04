
.PHONY: all

all: sgd_disc sgd_cont #pcor

# what about *.h?
#
CCFLAGS = -Wall -g -pg -std=gnu99
LIBRARIES = -lm

sgd_cont: common.c gmatrix.c scale.c loss.c evaluation.c sgd.c sgd_main.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o sgd_cont

sgd_disc: common.c gmatrix.c scale.c loss.c evaluation.c sgd.c sgd_main.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o sgd_disc -DDISCRETE

pcor: common.c gmatrix.c scale.c loss.c evaluation.c sgd.c pcor.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f sgd_cont sgd_disc pcor
# DO NOT DELETE
