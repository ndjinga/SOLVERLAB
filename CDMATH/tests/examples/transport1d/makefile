# Transport1d makefile

CC = g++
# Adapt the following line to your own system:
CDMATHDIR = ../../../../..
IFLAG = -I$(CDMATHDIR)/include -I.
LFLAG = -L$(CDMATHDIR)/lib
LIBS  =-linterpkernel -lmedC -lmedloader -lmedcoupling -lbase -lmesh -llinearsolver
OBJ = main.o

all: $(OBJ)
	$(CC) -o main $^ $(IFLAG) $(LFLAG) $(LIBS)
	
%.o: %.cxx
	$(CC) -c -o $@ $< $(CFLAGS) $(IFLAG) $(LFLAG) $(LIBS)
	
.PHONY: clean

clean:
	rm -f *.o *~ core $(INCDIR)/*~

sweep:
	rm -f *.vtu
	rm -f *.pvd
	rm -f *.csv
	rm -f *.med
