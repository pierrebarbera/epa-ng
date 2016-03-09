BIN=../bin
TEST=../test/data
LOGDIR=../test/log
OUTDIR=/tmp/epa
F+=
CF+=
TESTINPUT= $(TEST)/ref.tre $(TEST)/range_combined.fasta -g -w $(OUTDIR) $(F)

CC=clang++
#CC=g++
DEFS=
WARNINGS= -Wall -Wextra -pedantic -pedantic-errors
DBG=-g
CFLAGS= $(WARNINGS) -std=c++11 $(DBG) $(CF)
LINKFLAGS=$(DBG) #$(CFLAGS)

TARGET=$(BIN)/epamk

LIBS=-lm -lpll -lpll_optimize -lpll_msa -lpll_tree
CPP=$(wildcard *.cpp)
OBJS=$(notdir $(CPP:.cpp=.o))
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

callgrind: CFLAGS= $(WARNINGS) -std=c++11 -O2
profile: CFLAGS= $(WARNINGS) -std=c++11 -O2
# omp: CC=g++
omp: DBG=
omp: LIBS+= -liomp5
omp: CFLAGS+=-fopenmp=libomp
ranged_test: TESTINPUT= $(TEST)/ref.tre $(TEST)/range_combined.fasta -w $(OUTDIR) $(F)

all:  $(TARGET)

omp:  $(TARGET)

test: $(TARGET)
	mkdir -p $(OUTDIR)
	rm $(OUTDIR)/*
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

.PHONY: clean

clean:
	rm -f $(BIN)/*.o
	rm -f $(TARGET)
