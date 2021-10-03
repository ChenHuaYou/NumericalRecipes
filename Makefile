IDIR =include
CC=gcc
CFLAGS=-I$(IDIR)


_DEPS = nr.h nrutil.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gauss: xgaussj.o gaussj.o nrutil.o
	$(CC) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o *~ core
