SHELL = /bin/bash
CC ?= cc
FLAGS = -O3 -Wall -std=gnu11
LDFLAGS =

ARC = ../wdsat_solver

SRC_GEN = main.c wdsat_utils.c dimacs.c cnf.c xorset.c xorgauss.c wdsat.c fp.c
OBJ_GEN = $(SRC_GEN:.c=.o)
HDRS = config.h wdsat_utils.h dimacs.h cnf.h xorset.h xorgauss.h wdsat.h fp.h

### --- All Binaries

all: $(ARC)

s3-chain-19:
	$(MAKE) clean
	$(MAKE) ARC=../wdsat_s3_chain_19 FLAGS='$(FLAGS) -DWDSAT_PROFILE_S3_CHAIN_19'
	$(MAKE) clean

$(ARC): $(OBJ_GEN)
	$(CC) $^ -o $@ $(LDFLAGS) -lm

%.o: %.c $(HDRS)
	$(CC) $(FLAGS) -c $< -o $@

### --- OTHERS

clean:
	rm -f *.o *~ .*~

.PHONY: all clean s3-chain-19
