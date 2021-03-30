SOURCES = ballAlg.c ballAlg-omp.c gen_points.c ballQuery.c
OBJS = $(SOURCES:%.c=%.o)
CC = gcc

ifdef DEBUG
CFLAGS = -Wall -O3 -g -DDEBUG -fopenmp
else
CFLAGS = -Wall -O3 -fopenmp
endif

LDFLAGS = -lm
SERIAL = ballAlg
OMP = ballAlg-omp
TARGETS = $(SERIAL) $(OMP) ballQuery

all: $(TARGETS)


ballQuery: ballQuery.o
ballAlg: ballAlg.o gen_points.o
ballAlg-omp: ballAlg-omp.o gen_points.o

$(TARGETS):
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


ballQuery.o: ballQuery.c
ballAlg.o: ballAlg.c gen_points.h
gen_points.o: gen_points.c
ballAlg-omp.o: ballAlg-omp.c gen_points.h

$(OBJS):
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	@echo Cleaning...
	rm -f $(OBJS) $(TARGETS)

test: $(TARGETS)
	./test.sh
