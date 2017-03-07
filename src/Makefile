BIN=../bin
TEST=../test/data
LOGDIR=../test/log
OUTDIR=/tmp/epa
CEREAL=../libs/cereal/include/
CXXOPTS=../libs/cxxopts/src/
PLL_INC=../libs/include/libpll/
PLL_LIB=../libs/lib/

F+=
CF+=
BINARY_WRITE= -t $(TEST)/ref.tre -s $(TEST)/aln.fasta -q $(TEST)/query.fasta -B -w $(OUTDIR) $(F)
BINARY_READ=-b $(OUTDIR)/epa_binary_file -q $(TEST)/query.fasta -w $(OUTDIR) --filter-min-lwr 0.0 $(F)
NORM_TEST=-t $(TEST)/ref.tre -s $(TEST)/aln.fasta -q $(TEST)/query.fasta -w $(OUTDIR) $(F)

#CC=clang++
CC=g++
DEFS=
WARNINGS= -Wall -Wextra
DBG=-g

INCLUDES=-I$(CEREAL) -I$(CXXOPTS) -I$(PLL_INC)
LIB_PATHS=  -L$(PLL_LIB)
LIBS=$(PLL_LIB)libpll_optimize.a $(PLL_LIB)libpll_msa.a $(PLL_LIB)libpll_binary.a $(PLL_LIB)libpll.a -lm

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
omp: CFLAGS+=-fopenmp -D__OMP
ranged_test: BINARY_WRITE= -t $(TEST)/ref.tre -s $(TEST)/range_combined.fasta -r -w $(OUTDIR) $(F)
mpi: CC=skin mpicxx
mpi: CFLAGS+=-D__MPI
aa_test: BINARY_WRITE=-t $(TEST)/ref.tre -s $(TEST)/AA_aln.fasta -B -w $(OUTDIR) -m AA-GTR-DAYHOFF $(F)
aa_test: BINARY_READ=-b $(OUTDIR)/epa_binary_file -q $(TEST)/AA_query.fasta -w $(OUTDIR) -m AA-GTR-DAYHOFF --filter-min-lwr 0.0 $(F)
aa_test: NORM_TEST=-t $(TEST)/ref.tre -s $(TEST)/AA_aln.fasta -q $(TEST)/AA_query.fasta -w $(OUTDIR) -m AA-GTR-DAYHOFF $(F)

all:  $(TARGET)

omp:  $(TARGET)

mpi: $(TARGET)

test: $(TARGET)
	mkdir -p $(OUTDIR)
	-rm -f $(OUTDIR)/*
	$(TARGET) $(BINARY_WRITE)
	$(TARGET) $(BINARY_READ)
	$(TARGET) $(NORM_TEST)

aa_test: $(TARGET)
	mkdir -p $(OUTDIR)
	-rm -f $(OUTDIR)/*
	$(TARGET) $(BINARY_WRITE)
	$(TARGET) $(BINARY_READ)
	$(TARGET) $(NORM_TEST)

ranged_test: $(TARGET)
	mkdir -p $(OUTDIR)
	-rm $(OUTDIR)/*
	$(TARGET) $(BINARY_WRITE)
	$(TARGET) $(BINARY_READ)

mpi_test_valgrind: $(mpi)
	mkdir -p $(OUTDIR)
	-rm $(OUTDIR)/*
	-rm $(LOGDIR)/mpi.valgrind.*
	valgrind $(TARGET) $(BINARY_WRITE)
	mpirun -np 4 valgrind --log-file=$(LOGDIR)/mpi.valgrind.%p $(TARGET) $(BINARY_READ)

mpi_test: $(mpi)
	mkdir -p $(OUTDIR)
	-rm $(OUTDIR)/*
	$(TARGET) $(BINARY_WRITE)
	mpiexec -n 4 $(TARGET) $(BINARY_READ)


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
