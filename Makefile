SOURCES = ballAlg.c gen_points.c
OBJS = $(SOURCES:%.c=%.o)
CC = gcc
CFLAGS = -Wall
LDFLAGS = -fopenmp
SERIAL = ballAlg
TARGETS = $(SERIAL) ballQuery

all: $(TARGETS)


ballAlg: ballAlg.o gen_points.o
	$(CC) $(CFLAGS) $^ -o $@ -O3 $(LDFLAGS)


ballAlg.o: ballAlg.c gen_points.h
gen_points.o: gen_points.c

$(OBJS):
	$(CC) $(CFLAGS) -c $< -o $@


ballQuery: ballQuery.c
	$(CC) $(CFLAGS) $^ -o $@ -lm

clean:
	@echo Cleaning...
	rm $(OBJS) $(TARGETS)
