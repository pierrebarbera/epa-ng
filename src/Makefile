BIN=../bin
TEST=../test/data
LOGDIR=../test/log
OUTDIR=/tmp/epa
F+=
TESTINPUT= $(TEST)/ref.tre $(TEST)/aln.fasta -q $(TEST)/query.fasta -w $(OUTDIR) $(F)

CC=clang++
DEFS=
WARNINGS= -Wall -Wextra -pedantic -pedantic-errors
DBG=-g
CFLAGS= $(WARNINGS) -std=c++11 $(DBG)
LINKFLAGS=

TARGET=$(BIN)/epamk

LIBS=-lm -lpll -lpll_optimize
CPP=$(wildcard *.cpp)
OBJS=$(notdir $(CPP:.cpp=.o))
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

callgrind: CFLAGS= $(WARNINGS) -std=c++11 -O2

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
