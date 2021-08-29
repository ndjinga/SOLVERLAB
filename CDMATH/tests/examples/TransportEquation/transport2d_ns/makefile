# Transport2D_NS makefile

CC = g++
# Adapt the following line to your own system:
CDMATHDIR = ../../../../..
IFLAG = -I$(CDMATHDIR)/include -I.
LFLAG = -L$(CDMATHDIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib
LIBS  = -linterpkernel -lmedC -lmedloader -lmedcoupling -lbase -lmesh -llinearsolver -lpetsc

all: main main2

main: main.cxx 
	$(CC) -g -o $@ $^ $(IFLAG) $(LFLAG) $(LIBS)

main2: main2.cxx 
	$(CC) -g -o $@ $^ $(IFLAG) $(LFLAG) $(LIBS)

%.o: %.cxx
	$(CC) -c -g -o $@ $< $(CFLAGS) $(IFLAG) $(LFLAG) $(LIBS)

.PHONY: clean

clean:
	rm -f main main2 *.o *~ core $(INCDIR)/*~

sweep:
	rm -f *.vtu
	rm -f *.pvd
	rm -f *.csv
	rm -f Exercie2*.med
