
.PHONY: all

all: cd scale cd_double #pcor

# what about *.h?
#
#CCFLAGS = -Wall -g3 -std=gnu99
CCFLAGS = -Wall -g -pg -std=gnu99
LIBRARIES = -lm

cd: common.c gmatrix.c loss.c evaluation.c util.c cd.c main.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o cd

cd_double: common.c gmatrix.c loss.c evaluation.c util.c cd.c main.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o cd_double

scale: common.c gmatrix.c loss.c evaluation.c util.c cd.c scale.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o scale

pcor: common.c gmatrix.c loss.c evaluation.c util.c cd.c pcor.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f $(all)
