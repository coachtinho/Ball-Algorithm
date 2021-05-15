SOURCES = ballAlg.c ballAlg-omp.c ballAlg-mpi.c  gen_points.c ballQuery.c
OBJS = $(SOURCES:%.c=%.o)
CC = gcc
MPIC = mpicc

ifdef DEBUG
CFLAGS = -Wall -O3 -g -DDEBUG -fopenmp
else
CFLAGS = -Wall -O3 -g -fopenmp
endif

LDFLAGS = -lm
SERIAL = ballAlg
OMP = ballAlg-omp
MPI = ballAlg-mpi
TARGETS = $(SERIAL) $(OMP) $(MPI) ballQuery

all: $(TARGETS)

ballQuery: ballQuery.o
ballAlg: ballAlg.o gen_points.o
ballAlg-omp: ballAlg-omp.o gen_points.o
ballAlg-mpi: ballAlg-mpi.o gen_points.o

ballQuery:
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(SERIAL):
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(OMP):
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(MPI):
	$(MPIC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

ballQuery.o: ballQuery.c
ballAlg.o: ballAlg.c gen_points.h
gen_points.o: gen_points.c
ballAlg-omp.o: ballAlg-omp.c gen_points.h
ballAlg-mpi.o: ballAlg-mpi.c gen_points.h

$(OBJS):
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	@echo Cleaning...
	rm -f $(OBJS) $(TARGETS)

test: $(TARGETS)
	./test.sh
