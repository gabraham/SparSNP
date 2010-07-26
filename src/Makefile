
.PHONY: all


targets = cd scale cd_double
all: $(targets)

# what about *.h?
#
CCFLAGS = -Wall -pg -g3 -std=gnu99
#CCFLAGS = -Wall -O2 -std=gnu99
LIBRARIES = -lm

cd: common.c gmatrix.c loss.c util.c main.c cd.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o cd

cd_double: common.c gmatrix.c loss.c util.c main.c cd.c
	gcc -DDTYPE=double $(CCFLAGS) $(LIBRARIES) $^ -o cd_double

scale: common.c gmatrix.c loss.c util.c cd.c scale.c
	gcc $(CCFLAGS) $(LIBRARIES) $^ -o scale

pcor: common.c gmatrix.c loss.c util.c cd.c pcor.c
	gcc -Dintype=double $(CCFLAGS) $(LIBRARIES) $^ -o pcor

clean:
	/bin/rm -f $(targets)
