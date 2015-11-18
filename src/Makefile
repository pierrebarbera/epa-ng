BIN=../bin
TEST=../test/data
TESTINPUT= $(TEST)/ref.tre $(TEST)/aln.fasta $(TEST)/query.fasta

CC=clang++
DEFS=
WARNINGS= -Wall -Wextra -pedantic -pedantic-errors
CFLAGS= $(WARNINGS) -std=c++11 -g
LINKFLAGS=

TARGET=$(BIN)/epa

LIBS=-lm -lpll
OBJS=Main.o MSA.o Tree.o Model.o file_io.o Sequence.o
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
