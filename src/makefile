CC=g++
CFLAGS=-std=c++11 -g -Wall
OBJECTS=main.o basisState.o projector.o hamiltonian.o
NAME=fciqmctest

all: $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS) $(CFLAGS)

$(OBJECTS): %.o: %.cpp $(TEMPLATES)
	$(CC) -c $(CFLAGS) $< -o $@
