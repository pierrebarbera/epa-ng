BIN=../bin

CC=clang++
DEFS=
WARNINGS= -Wall -Wextra -Wno-c++98-compat
CFLAGS= $(WARNINGS) -std=c++11
LINKFLAGS= -L/usr/local/lib/

TARGET=$(BIN)/epa

LIBS=-lpll -lm
OBJS=Main.o MSA.o Tree.o 
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

all:  $(TARGET)

run: $(TARGET)
	$(TARGET)

$(TARGET): $(OUT_OBJS)
	$(CC) $(LINKFLAGS) $(LIBS) -o $@ $+ 

$(OUT_OBJS): $(BIN)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean 

clean:
	rm -f $(BIN)/*.o
	rm -f $(TARGET)