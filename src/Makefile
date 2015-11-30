BIN=../bin
TEST=../test/data
TESTINPUT= $(TEST)/ref.tre $(TEST)/aln.fasta -q $(TEST)/query.fasta

CC=clang++
DEFS=
WARNINGS= -Wall -Wextra -pedantic -pedantic-errors
CFLAGS= $(WARNINGS) -std=c++11 -g
LINKFLAGS=

TARGET=$(BIN)/epamk

LIBS=-lm -lpll -lpll_optimize
CPP=$(wildcard *.cpp)
OBJS=$(notdir $(CPP:.cpp=.o))
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

all:  $(TARGET)

test: $(TARGET)
	$(TARGET) $(TESTINPUT)

valgrind: $(TARGET)
	valgrind $(TARGET) $(TESTINPUT)

gdb: $(TARGET)
	gdb --args $(TARGET) $(TESTINPUT)

$(TARGET): $(OUT_OBJS)
	$(CC) $(LINKFLAGS) -o $@ $+ $(LIBS)

$(OUT_OBJS): $(BIN)/%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean

clean:
	rm -f $(BIN)/*.o
	rm -f $(TARGET)
