FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: spmat.o main.o SetOfGroups.o modularity.o BHatMatrix.o
	gcc main.o spmat.o SetOfGroups.o modularity.o BHatMatrix.o -o cluster $(LIBS)
clean:
	rm -rf *.o cluster

spmat.o: spmat.c spmat.h
	gcc $(FLAGS) -c spmat.c

SetOfGroups.o: SetOfGroups.c SetOfGroups.h BHatMatrix.h
	 gcc $(FLAGS) -c SetOfGroups.c
 
BHatMatrix.o: BHatMatrix.c BHatMatrix.h spmat.h 
	gcc $(FLAGS) -c BHatMatrix.c

modularity.o: modularity.c modularity.h BHatMatrix.h SetOfGroups.h
	gcc $(FLAGS) -c modularity.c

main.o: main.c BHatMatrix.h SetOfGroups.h modularity.h
	gcc $(FLAGS) -c main.c
