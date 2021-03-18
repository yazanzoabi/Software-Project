#include <stdio.h>
#include <stdlib.h>
#include "SPBufferset.h"
#include "SetOfGroups.h"
#include "BHatMatrix.h"
#define EPSILON 0.00001

/*
 * struct _Vectors
 *
 * 	saves vectors throughout the use of functions,
 *  to eliminate the need to allocate more space than needed.
 */
typedef struct _Vectors {
	double **vectors;
	int* index;
	void (*free)(struct _Vectors* Vec);
} Vector;

/**
 * void free_Vec(struct _Vectors* Vec)
 *
 *  Frees all resources used by Vec.
 */
void free_Vec(struct _Vectors* Vec){
	free(Vec->vectors[0]);
	free(Vec->vectors[1]);
	free(Vec->vectors[2]);
	free(Vec->vectors);
	free(Vec->index);
	free(Vec);

}

/**
 * Vector* allocate_Vector(int n)
 *
 *  Allocates the initial space for 3 vectors of type double and one of type int with size n for each one,
 *  these three vectors will be used throughout the project,
 *  the first vector represents s.
 *  the other vectors will be used for different other uses.
 *
 */
Vector* allocate_Vector(int n){
	Vector* Vec=malloc(sizeof(Vector));
	if(Vec ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	Vec->vectors=malloc(sizeof(double*)*3);
	if(Vec->vectors ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	Vec->vectors[0]=malloc(sizeof(double)*n);
	if(n != 0 && Vec->vectors[0] ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	Vec->vectors[1]=malloc(sizeof(double)*n);
	if(n != 0 && Vec->vectors[1] ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	Vec->vectors[2]=malloc(sizeof(double)*n);
	if(n != 0 && Vec->vectors[2] ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	Vec->index=malloc(sizeof(int)*n);
	if(n != 0 && Vec->index ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	Vec->free=free_Vec;

	return Vec;
}

/**
 * void Random_vector(double *vector, int n)
 *
 * generates a random vector of doubles at range (0,10) of size n s.t. vector[0] != 0.
 * the generated vector is written on the RandVec.
 * (RandVec is pre-allocated)
 */
void Random_vector(double *RandVec, int n){
	register double *ptr1;
	register int i;

	ptr1 = RandVec;
	for(i = 0; i < n; i++){
		if(i != 01)
			*ptr1 = rand()%10;
		else{
			*ptr1 = rand()%10;
			while(*ptr1 == 0)
				*ptr1 = rand()%10;
		}
		ptr1++;
	}
}

/*
 * int max(double *a,int n)
 *
 * returns the index of max value in a s.t. a with length n
 * , and n-1 if no positive value is found.
 */
int max(double *a,int n){
	register double max = *a;
	register int i, index = 0;

	for (i = 0; i < n; i++){
		if(max <= *a){
			index = i;
			max = *a;
		}
		a++;
	}
	if (max <= 0){
		return n-1;
	}
	return index;
}

/*
 * void buildSet(groups* Set, double *s,group *g)
 *
 * group g, s vector include 1 -1 s.t. represent division of g,
 * allocate a new group according to a division into Set.
 *(Set is pre-allocated)
 */
void buildSet(groups* Set, double *s,group *g){
	group *g1, *g2;
	double *sPtr;
	int *indexes1Ptr, *indexes2Ptr,*indexes;
	int i;
	int* indexes1;
	int* indexes2;
	int size1 = 0, size2 = 0;
	BHat** Bhats = NULL;
	int size;


	size=g->size;
	sPtr = s;

	/*compute sizes of groups*/
	for(i = 0; i < size; i++){
		if(*sPtr == 1)
			size1++;
		if(*sPtr == -1)
			size2++;
		sPtr++;
	}


	if(size1 == 0 || size2 == 0){
		Set->Insert(Set, g);
		return;
	}

	indexes=g->indexes;
	/*compute indexes of each group*/
	indexes1=malloc(sizeof(int) * size1);
	if(size1 != 0 && indexes1 ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	indexes2=malloc(sizeof(int) * size2);
	if(size2 != 0 && indexes2 ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	sPtr = s;
	indexes1Ptr = indexes1;
	indexes2Ptr = indexes2;
	for(i = 0; i < size; i++){
			if(*sPtr == 1){
				*indexes1Ptr = *indexes;
				indexes1Ptr++;
			}
			if(*sPtr == -1){
				*indexes2Ptr = *indexes;
				indexes2Ptr++;
			}
			sPtr++;
			indexes++;
	}

	/*Build division of g*/
	Bhats = Compute_subBHat( g->B, s, size1, size2);
	g1 = group_allocate(size1, indexes1, Bhats[0]);
	g2 = group_allocate(size2, indexes2, Bhats[1]);
	Set->Insert(Set, g1);
	Set->Insert(Set, g2);

	free(Bhats);
	g->free(g);
}

/*
 * void maximization(Vector *vec, group* g)
 *
 *  Optimizes the division of g, by using smart shortcuts and optimizations.
 *
 *  1. a vec pointer with the s in vec.vectors[0] calculated by algorithm 2,
 *  2. and a group pointer with the relevant B matrix.
 *
 *  the new division is then written to s vec.
 *
 */
void maximization(Vector *vec, group* g){
	double *score, *s, *x, *improve;
	double delta = EPSILON;
	double *IPtr, *scorePtr, *sPtr, *xPtr;
	double maxDelta;
	double **vectores;
	register int *improveIndex, *moved, *degrees;
	register int i, j, first, M;
	int maxIndex;
	int *IIPtr, *movedPtr, *degreesPtr;
	int n = g->size;
	BHat* B;

	B = g->B;
	M = B->M;

	vectores = vec->vectors;
	s = *vectores;
	score = *(vectores +1);
	improve = *(vectores +2);
	improveIndex = vec->index;
	degrees = B->degrees;
	x = malloc(sizeof(double)*n);
	if(n != 0 && x ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	/* initial score values*/
	B->B_mult(B, s, x);
	sPtr = s; xPtr = x; degreesPtr = degrees; scorePtr = score;
	for(i=0; i < n; i++){
		*scorePtr = -2 * ((*sPtr * *xPtr) + ((double)*degreesPtr * *degreesPtr / M));
		sPtr++; xPtr++; degreesPtr++; scorePtr++;
	}

	/*REPEAT*/
	while (delta >= EPSILON){

		moved = (int*)calloc(n, sizeof(int));
		if(n != 0 && moved == NULL){
			printf("A fatal memory error has been detected.. exiting.\n");
			exit(-1);
		}
		movedPtr = moved;

		IPtr = improve;
		IIPtr = improveIndex;

		movedPtr = moved;
		for(j = 0; j < n; j++){
			movedPtr = moved;
			first = 0;
			scorePtr = score;
			/*find max delta*/
			for(i = 0; i < n; i++){
				if(*movedPtr == 0){
					delta = *scorePtr;
					if(first == 0){
						first = 1;
						maxDelta = delta;
						maxIndex = i;
					}
					else if(maxDelta < delta){
						maxDelta = delta;
						maxIndex = i;
					}
				}
				scorePtr++;
				movedPtr++;
			}

			moved[maxIndex] = 1;
			s[maxIndex] = -s[maxIndex];
			B->BupdateModularity(B, score, s, s[maxIndex], maxIndex);

			*IIPtr = maxIndex;
			IIPtr++;

			if(j == 0){
				*IPtr = maxDelta;
				IPtr++;
			}
			else{
				*IPtr = maxDelta + *(IPtr -1);
				IPtr++;
			}

		}

		maxIndex = max(improve,n);
		IIPtr--;

		for(i = n-1; i > maxIndex; i--){
			j = *IIPtr;
			s[j] = -s[j];
			B->BupdateModularity(B, score, s, s[j], j);
			IIPtr--;
		}

		if(maxIndex == n-1)
			delta = 0;
		else
			delta = improve[maxIndex];

		free(moved);
	}
	free(x);
}

/*
 * void algorithm2( group *g,Vector* vec)
 *
 * Devises the group g by the eigen vector.
 * and writes the division into s in vec->vectors[0],
 */
void algorithm2( group *g, Vector* vec){
	register double *sPtr, *resultPtr;
	register double eig, mod;
	register int sizeofgroup, i;
	BHat* Bmatrix;
	double** vectores;
	double *RandVec, *s, *result;


	sizeofgroup = g->size;

	Bmatrix = g->B;
	vectores = vec->vectors;
	s = *vectores;
	RandVec = *(vectores +1);
	result = *(vectores +2);
	Random_vector(RandVec, sizeofgroup);

	/*compute the eigenvector and eigenvalue*/
 	eig = Bmatrix->get_eigenValue(Bmatrix, RandVec, result);

	if(eig <= EPSILON){
		sPtr = s;
		for(i = 0; i < sizeofgroup; i++){
			*sPtr = 1;
			sPtr++;
		}

		maximization(vec, g);
		return;
	}

	/* builds the initial division in s.*/
	sPtr = s;
	resultPtr = result;
	for(i = 0; i < sizeofgroup; i++){

		if(*resultPtr > EPSILON)
			*sPtr = 1;
		else
			*sPtr = -1;
		sPtr++;
		resultPtr++;
	}


	/*compute sBHat[g]s */
	Bmatrix->BHat_mult(Bmatrix, s, result);
	mod = Vector_mult(s, result, sizeofgroup);

	if(mod <= EPSILON){
		sPtr = s;
		for(i = 0; i < sizeofgroup; i++){
			*sPtr = 1;
			sPtr++;
		}
	}

}


/*
 * groups* algorithm3(BHat *B)
 *
 * returns Set of groups with the complete division of the BHat matrix B.
 */

groups* algorithm3(BHat *B){

	register int i, n;
	register groups *P,*O,*Set = NULL;
	group *g1, *g2 = NULL, *g3 = NULL, *g4 = NULL;
	Vector* Vec;
	int* indexes;

	P = SetOfGroups_allocate();
	O = SetOfGroups_allocate();
	Set = SetOfGroups_allocate();
	n = B->n;

	Vec = allocate_Vector(n);
	indexes = malloc(sizeof(int)*n);
	if(n != 0 && indexes ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}


	/*trivial division*/
	for(i = 0; i< n ; i++){
		indexes[i] = i;
	}

	g1 = group_allocate(n, indexes, B);

	P->Insert(P, g1);

	while(P->size != 0){
		g2 = P->Remove(P);

		algorithm2(g2, Vec);
		maximization(Vec, g2);
		buildSet(Set, *(Vec->vectors), g2);

		if(Set->size == 1){
			g2 = Set->Remove(Set);
			O->Insert(O, g2);
		}
		else{
			g3 = Set->Remove(Set);
			g4 = Set->Remove(Set);

			if(g3->size == 1)
				O->Insert(O, g3);
			else
				P->Insert(P, g3);

			if(g4->size == 1)
				O->Insert(O, g4);
			else
				P->Insert(P, g4);

		}
	}

	P->free(P);
	Set->free(Set);
	Vec->free(Vec);
	return O;
}
