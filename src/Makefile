BIN=../bin
TEST=../test/data
LOGDIR=../test/log
OUTDIR=/tmp/epa
CEREAL=../libs/cereal/include/
PLL=../libs/pll-modules/libs/libpll/src/
PLL_BINARY=../libs/pll-modules/src/binary/
PLL_OPTIMIZE=../libs/pll-modules/src/optimize/
PLL_MSA=../libs/pll-modules/src/msa/

F+=
CF+=
BINARY_WRITE= -t $(TEST)/ref.tre -s $(TEST)/aln.fasta -q $(TEST)/query.fasta -g -B -w $(OUTDIR) $(F)
BINARY_READ=-b $(OUTDIR)/epa_binary_file -q $(TEST)/query.fasta -w $(OUTDIR)
NORM_TEST=-t $(TEST)/ref.tre -s $(TEST)/aln.fasta -q $(TEST)/query.fasta -g -w $(OUTDIR) $(F)

#CC=clang++
CC=g++
DEFS=
WARNINGS= -Wall -Wextra
DBG=-g

INCLUDES=-I$(CEREAL) -I$(PLL) -I$(PLL_BINARY) -I$(PLL_OPTIMIZE) -I$(PLL_MSA)
RPATHS= -Wl,-rpath,$(PLL).libs/,-rpath,$(PLL_BINARY),-rpath,$(PLL_OPTIMIZE),-rpath,$(PLL_MSA)
LIB_PATHS=  -L$(PLL).libs/ -L$(PLL_BINARY) -L$(PLL_OPTIMIZE) -L$(PLL_MSA)
LIBS=-lm -lpll -lpll_optimize -lpll_msa -lpll_binary

CFLAGS= $(WARNINGS) -std=c++11 $(DBG) $(INCLUDES) $(CF)
LINKFLAGS=$(DBG) $(RPATHS) $(LIB_PATHS)

TARGET=$(BIN)/epamk

CPP=$(wildcard *.cpp)
OBJS=$(notdir $(CPP:.cpp=.o))
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

callgrind: CFLAGS= $(WARNINGS) -std=c++11 -O2
profile: CFLAGS= $(WARNINGS) -std=c++11 -O2
omp: CC=g++
omp: DBG=
omp: LIBS+=-lgomp
omp: CFLAGS+=-fopenmp
ranged_test: BINARY_WRITE= -t $(TEST)/ref.tre -s $(TEST)/range_combined.fasta -r -w $(OUTDIR) $(F)
mpi: CC=mpic++
mpi: CFLAGS+=-D__MPI

all:  $(TARGET)

omp:  $(TARGET)

mpi: $(TARGET)

test: $(TARGET)
	mkdir -p $(OUTDIR)
	-rm -f $(OUTDIR)/*
	$(TARGET) $(BINARY_WRITE)
	$(TARGET) $(BINARY_READ)

ranged_test: $(TARGET)
	mkdir -p $(OUTDIR)
	-rm $(OUTDIR)/*
	$(TARGET) $(BINARY_WRITE)
	$(TARGET) $(BINARY_READ)

mpi_test_valgrind: $(mpi)
	mkdir -p $(OUTDIR)
	-rm $(OUTDIR)/*
	-rm $(LOGDIR)/mpi.valgrind.*
	$(TARGET) $(BINARY_WRITE)
	mpirun -np 4 valgrind --log-file=$(LOGDIR)/mpi.valgrind.%p $(TARGET) $(BINARY_READ)

mpi_test: $(mpi)
	mkdir -p $(OUTDIR)
	-rm $(OUTDIR)/*
	$(TARGET) $(BINARY_WRITE)
	mpirun -np 4 $(TARGET) $(BINARY_READ)


valgrind: $(TARGET)
	mkdir -p $(OUTDIR)
	valgrind $(TARGET) $(NORM_TEST)
	valgrind $(TARGET) $(BINARY_WRITE)
	valgrind $(TARGET) $(BINARY_READ)

callgrind: $(TARGET)
	mkdir -p $(LOGDIR)
	mkdir -p $(OUTDIR)
	valgrind --tool=callgrind --callgrind-out-file=$(LOGDIR)/callgrind.output.%p \
	$(TARGET) $(BINARY_WRITE)

gdb: $(TARGET)
	mkdir -p $(OUTDIR)
	gdb --args $(TARGET) $(NORM_TEST)

$(TARGET): $(OUT_OBJS)
	$(CC) $(LINKFLAGS) -o $@ $+ $(LIBS)

$(OUT_OBJS): $(BIN)/%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean, mpi, omp, all, test, ranged_test, valgrind, callgrind, gdb

clean:
	-rm -f $(BIN)/*.o
	-rm -f $(TARGET)
