BIN=../bin

CC=clang++
DEFS=
WARNINGS= -Weverything -Wno-c++98-compat
CFLAGS= $(WARNINGS) -I../../libpll/src -std=c++11
LINKFLAGS= -L/usr/local/lib/

TARGET=$(BIN)/epa

LIBS=-lpll -lm
OBJS=epa.o Main.o 
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