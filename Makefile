OBJ    = main.o solverwa.o solver.o mq.o io.o mff.o rhyd.o eos.o matrix.o
CC=g++
CFLAGS=-W -Wall -Wextra -std=c++11 -pedantic -O3 -c
all: exe

exe: $(OBJ)
	$(CC) -fopenmp -o exe $(OBJ)

main.o: main.cpp solver.hpp solverwa.hpp mq.hpp
	$(CC) -fopenmp $(CFLAGS) main.cpp

solver.o: solver.cpp solver.hpp mff.hpp io.hpp
	$(CC) $(CFLAGS) solver.cpp

solverwa.o: solverwa.cpp solverwa.hpp solver.hpp
	$(CC) $(CFLAGS) solverwa.cpp

mq.o: mq.cpp mq.hpp matrix.hpp io.hpp
	$(CC) $(CFLAGS) mq.cpp

io.o: io.cpp io.hpp matrix.hpp
	$(CC) $(CFLAGS) io.cpp

mff.o: mff.cpp mff.hpp rhyd.hpp
	$(CC) $(CFLAGS) mff.cpp

rhyd.o: rhyd.cpp matrix.hpp rhyd.hpp eos.hpp
	$(CC) $(CFLAGS) rhyd.cpp

eos.o: eos.cpp eos.hpp matrix.hpp
	$(CC) $(CFLAGS) eos.cpp

matrix.o: matrix.cpp matrix.hpp
	$(CC) $(CFLAGS) matrix.cpp

