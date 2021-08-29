# thermique1d makefile

CC = g++

# Adapt the following lines to your own system:
PETSCDIR = /usr/lib/petscdir/3.4.5
MPIDIR = /usr/lib/openmpi
CDMATHDIR = ../../../../..

IFLAG = -I$(PETSCDIR)/include -I$(MPIDIR)/include -I$(CDMATHDIR)/include -I.
LFLAG = -L$(CDMATHDIR)/lib
LIBS  = -linterpkernel -lmedC -lmedloader -lmedcoupling -lbase -lmesh -llinearsolver
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
