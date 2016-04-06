BIN=../bin
TEST=../test/data
LOGDIR=../test/log
OUTDIR=/tmp/epa
CEREAL=../../cereal/include/
F+=
CF+=
TESTINPUT= -t $(TEST)/ref.tre -s $(TEST)/aln.fasta -q $(TEST)/query.fasta -g -w $(OUTDIR) $(F)

CC=clang++
#CC=g++
DEFS=
WARNINGS= -Wall -Wextra -pedantic -pedantic-errors
DBG=-g
CFLAGS= $(WARNINGS) -std=c++11 $(DBG) -I$(CEREAL) $(CF)
LINKFLAGS=$(DBG) #$(CFLAGS)

TARGET=$(BIN)/epamk

LIBS=-lm -lpll -lpll_optimize -lpll_msa -lpll_tree -lpll_binary
CPP=$(wildcard *.cpp)
OBJS=$(notdir $(CPP:.cpp=.o))
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

callgrind: CFLAGS= $(WARNINGS) -std=c++11 -O2
profile: CFLAGS= $(WARNINGS) -std=c++11 -O2
# omp: CC=g++
omp: DBG=
omp: LIBS+= -liomp5
omp: CFLAGS+=-fopenmp=libomp
ranged_test: TESTINPUT= -t $(TEST)/ref.tre -s $(TEST)/range_combined.fasta -r -w $(OUTDIR) $(F)
mpi: CC=mpic++
mpi: CFLAGS+=-D__MPI

all:  $(TARGET)

omp:  $(TARGET)

mpi: $(TARGET)

test: $(TARGET)
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*
	$(TARGET) $(TESTINPUT)

ranged_test: $(TARGET)
	mkdir -p $(OUTDIR)
	rm $(OUTDIR)/*
	$(TARGET) $(TESTINPUT)


valgrind: $(TARGET)
	mkdir -p $(OUTDIR)
	valgrind $(TARGET) $(TESTINPUT)

callgrind: $(TARGET)
	mkdir -p $(LOGDIR)
	mkdir -p $(OUTDIR)
	valgrind --tool=callgrind --callgrind-out-file=$(LOGDIR)/callgrind.output.%p \
	$(TARGET) $(TESTINPUT)

profile: $(TARGET)
	mkdir -p $(OUTDIR)
	$(TARGET) $(TESTINPUT) -b 10

gdb: $(TARGET)
	mkdir -p $(OUTDIR)
	gdb --args $(TARGET) $(TESTINPUT)

$(TARGET): $(OUT_OBJS)
	$(CC) $(LINKFLAGS) -o $@ $+ $(LIBS)

$(OUT_OBJS): $(BIN)/%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean, mpi, omp, all, test, ranged_test, valgrind, callgrind, gdb

clean:
	rm -f $(BIN)/*.o
	rm -f $(TARGET)
