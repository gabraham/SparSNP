
.PHONY: all

all: sgd_cont sgd_disc

# what about *.h?

sgd_cont: sgd.c gmatrix.c loss.c evaluation.c common.c
	gcc -Wall -O2 -g $^ -o sgd_cont

sgd_disc: sgd.c gmatrix.c loss.c evaluation.c common.c
	gcc -Wall -O2 -g $^ -o sgd_disc -DDISCRETE

clean:
	/bin/rm sgd_cont sgd_disc
