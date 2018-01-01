#Makefile for 2Dcluster

MAINSRC = main
OUT = a.out
OUTDIR = ./
ALLOBJ = main.o func.o
LIBS = -lm
LIBDIR = /sw/lib
CXX = gcc
CXXFLAGS = -O3 #-Wunused-variable -Wunused-but-set-variable

$(OUT): $(MAINSRC).o func.o
	$(CXX) -O3 -L$(LIBDIR) $(ALLOBJ) $(LIBS) -Wall -o $(OUTDIR)$(OUT)
	echo "Linked "$(OUT)" !"
	rm -f *.o

$(MAINSRC).o: $(MAINSRC).c main.h 
	$(CXX) -c $(CXXFLAGS) -Wall $(MAINSRC).c
	echo "Compiled main"

func.o: func.c main.h Mt.h
	$(CXX) -c $(CXXFLAGS) -Wall func.c
	echo "Compiled func"

clean: 
	rm -f $(OUT) 
	rm -f $(ALLOBJ)
	rm -f *.o
