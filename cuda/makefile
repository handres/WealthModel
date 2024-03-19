CC=g++
CFLAGS=-std=c++11 -lgsl -lgslcblas -O3


all: simple_plot.o
	$(CC) *.o -o simple_plot $(CFLAGS)
simple_plot.o: simple_plot.cpp
	$(CC) simple_plot.cpp -c -o simple_plot.o $(CFLAGS)
