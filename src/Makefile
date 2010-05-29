
.PHONY: all

all: sgd_cont sgd_disc pcor

# what about *.h?
#
CCFLAGS = -Wall -g -std=gnu99
LIBRARIES = -lm

sgd_cont: sgd_main.c sgd.c gmatrix.c loss.c evaluation.c common.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o sgd_cont

sgd_disc: sgd_main.c sgd.c gmatrix.c loss.c evaluation.c common.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o sgd_disc -DDISCRETE

pcor: pcor.c sgd.c gmatrix.c loss.c evaluation.c common.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm sgd_cont sgd_disc
