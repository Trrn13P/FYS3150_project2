CPPflags= c++ -o
LIB = -larmadillo -llapack -lblas

all: compile execute

compile:
	${CPPflags} ./main.out main.cpp eigenvalues.cpp ${LIB}

execute:
	./main.out
