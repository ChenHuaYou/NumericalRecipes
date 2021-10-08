CC=g++
CFLAGS=-I./include -g -fsanitize=undefined,address,leak
LFLAGS=-fsanitize=undefined,address,leak

_DEPS = nr.h
DEPS = $(patsubst %,include/%,$(_DEPS))

_OBJ = gaussj.o xgaussj.o
OBJ = $(patsubst %,obj/%,$(_OBJ))

obj/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gauss: obj/gaussj.o obj/xgaussj.o
	$(CC) -o $@ $^ $(LFLAGS)

.PHONY: clean

clean:
	rm -f obj/*.o *~ core gauss
