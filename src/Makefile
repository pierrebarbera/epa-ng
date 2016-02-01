BIN=../bin
TEST=../test/data
LOGDIR=../test/log
OUTDIR=/tmp/epa
F+=
CF+=
TESTINPUT= $(TEST)/ref.tre $(TEST)/aln.fasta -q $(TEST)/query.fasta -w $(OUTDIR) -oO $(F)

CC=clang++
#CC=g++
DEFS=
WARNINGS= -Wall -Wextra -pedantic -pedantic-errors
DBG=-g
CFLAGS= $(WARNINGS) -std=c++11 $(DBG) -Ofast $(CF)
LINKFLAGS=$(DBG) #$(CFLAGS)

TARGET=$(BIN)/epamk

LIBS=-lm -lpll -lpll_optimize
CPP=$(wildcard *.cpp)
OBJS=$(notdir $(CPP:.cpp=.o))
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

callgrind: CFLAGS= $(WARNINGS) -std=c++11 -O2
profile: CFLAGS= $(WARNINGS) -std=c++11 -O2
gprof: LINKFLAGS+= --static -pg
gprof: CFLAGS+= -pg

all:  $(TARGET)

test: $(TARGET)
	mkdir -p $(OUTDIR)
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

gprof: $(TARGET)
	mkdir -p $(LOGDIR)
	mkdir -p $(OUTDIR)
	$(TARGET) $(TESTINPUT)
	gprof $(OUTDIR)/gmon.out


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
