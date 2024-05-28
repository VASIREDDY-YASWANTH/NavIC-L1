IDIR =./include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =./lib

LIBS=-lm

_DEPS = position_main.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_SRC = position_main.c
SRC = $(patsubst %,src/%,$(_SRC))

_OBJ = $(_SRC:.c=.o) 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

position_calc: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# .PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 