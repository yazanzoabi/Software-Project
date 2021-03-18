
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPBufferset.h"
#include "spmat.h"
/**
 * typedef struct row
 *
 * represents a row in matrix spmat,
 * s.t pointer row include the numbers of columns of nnz values,
 * num_of_nnz represents the number of nnz values in row.
 */
typedef struct row {
	int num_of_nnz;
	int* row;
} NodeOfRow;

/**
 * typedef struct list_spmat
 *
 * array of linked lists of NodesOfRows.
 * which represents the values of spmat.
 */
typedef struct list_spmat {
	NodeOfRow** lst;
} ListOfNode;

/*
 * double shift_mat(struct _spmat *A, int *degrees,int M)
 *
 * calculates the shift amount of a matrix B, s.t B[i][j] = A[i][j] -degrees[i]*degrees[j]/M.
 *  this function is called once, the result is then used for all sub B matrixes.
 */
double shift_mat(struct _spmat *A, int *degrees,int M){
	register int j,i;
	register double sum, max = 0;
	register int *degreesPtr1, *degreesPtr2;
	NodeOfRow** lst = ((ListOfNode*)A->private)->lst;
	int* entry;
	NodeOfRow** row;
	register int num_col, index;
	register int n, sizeofrow;

	n = A->n;

	row = lst;
	degreesPtr1 = degrees;
	for(j = 0; j < n; j++){
		sum = 0;
		entry = (*row)->row;
		sizeofrow = (*row)->num_of_nnz;
		i = 0;
		degreesPtr2 = degrees;

		while(sizeofrow != 0 && i < n){
			num_col = *entry;
			index = i;
			if(index < num_col){
				sum += fabs(-(double)*degreesPtr1**degreesPtr2/M);
				i++;
				degreesPtr2++;
			}
			else{
				sum += fabs(1-(double)*degreesPtr1**degreesPtr2/M);
				degreesPtr2++;
				entry++;
				i++;
				sizeofrow--;
			}

		}

		for(;i < n; i++){
			sum += fabs(-(double)*degreesPtr1**degreesPtr2/M);

			degreesPtr2++;
		}

		if(max < sum)
			max = sum;
		row++;
		degreesPtr1++;
	}

	return max;
}

/*
 * void add_row(struct _spmat *A, int *row, int i,int size)
 *
 * Adds row i to the matrix. Called before any other call,
 * exactly n times in order (i = 0 to n-1)
 */
void add_row(struct _spmat *A, int *row, int i,int size){
	NodeOfRow* Row = NULL;
	Row = (NodeOfRow*)malloc(sizeof(NodeOfRow));

	if(Row ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	Row->num_of_nnz = size;
	Row->row = row;
	((ListOfNode*)A->private)->lst[i] = Row;
}

/*
 * void mat_free(struct _spmat *A)
 *
 * Frees all resources used by A
 */
void mat_free(struct _spmat *A){

	NodeOfRow** lst = ((ListOfNode*)A->private)->lst;
	NodeOfRow* RowPtr;

	int n = A->n, i;

	for(i = 0; i < n; i++){
		RowPtr = *lst;
		if(RowPtr->row != NULL)
			free(RowPtr->row);
		free(RowPtr);
		lst++;
	}
	free(((ListOfNode*)A->private)->lst);
	free(A->private);
	free(A);
}

/*
 * 	void mult(const struct _spmat *A, const double *v, double *result)
 *
 *  Multiplies matrix A by vector v,
 *  into result (result is pre-allocated)
 */
void mult(const struct _spmat *A, const double *v, double *result){

	register double sum;
	register int j, i, col, prevcol = -1;
	register int* entry;
	register NodeOfRow** row = ((ListOfNode*)A->private)->lst;
	register int n, sizeofrow;
	const double *vPtr;
	n = A->n;
	j = 0;

	for (j = 0; j < n; j++){
		entry = (*row)->row;
		sizeofrow = (*row)->num_of_nnz;
		vPtr=v;
		sum = A->shift_amount * v[j];
		for(i = 0; i < sizeofrow; i++){
			col = *entry;
			if(prevcol == -1)
				vPtr += col;
			else
				vPtr += col - prevcol;
			sum += *vPtr;
			entry++;
			prevcol = col;
		}
		prevcol = -1;
		row++;
		*result = sum;
		result++;
	}

}

/*
 * void Sum_OF_Rows(const spmat *A,double *sumOFRows)
 *
 * compute sum of each row in A,
 * into sumOFRows (sumOFRows pre-allocated).
 *
 */
void Sum_OF_Rows(const spmat *A,double *sumOFRows){
	register NodeOfRow** lst = ((ListOfNode*)A->private)->lst;
	register int j;
	register double sum;
	register int n;

	n = A->n;
	for(j = 0; j < n; j++){
		sum = (*lst)->num_of_nnz;
		*sumOFRows = *sumOFRows + sum;
		sumOFRows++;
		lst++;
	}
}

/*
 * void A_update_Modularity(const struct _spmat *A, double *result, double *s, int sIndex, int index)
 *
 *  updates score s.t:
 *   score[i] =-4*A[i][index]*s[i]*sIndex for i != index.
 */
void A_update_Modularity(const struct _spmat *A, double *result, double *s, int sIndex, int index){
	NodeOfRow** lst = ((ListOfNode*)A->private)->lst;
	int* entry;
	int col, prevcol = -1, n;

	n = lst[index]->num_of_nnz;
	if(n == 0){
		return;
	}
	entry = lst[index]->row;
	while(n!=0){
		col = *entry;
		if (prevcol==-1)
			result+=col;
		else
			result += col - prevcol;
		if(col != index){
			*result -= sIndex*4*s[col];

		}
		entry++;
		n--;
		prevcol = col;
	}
}

/*
 * spmat* spmat_allocate(int n)
 *
 * Allocates a new spmat sparse matrix of size n.
 */
spmat* spmat_allocate(int n){

	spmat* mat = malloc(sizeof(spmat));
	ListOfNode* list = malloc(sizeof(ListOfNode));

	if(mat == NULL || list ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	list->lst = malloc(sizeof(NodeOfRow*) * n);
	if(n != 0 && list->lst ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	mat->n = n;
	mat->private = list;
	mat->add_row = add_row;
	mat->free = mat_free;
	mat->mult = mult;
	mat->shift_amount = 0;
	mat->shift_mat = shift_mat;
	mat->Sum_OF_Rows = Sum_OF_Rows;
	mat->AupdateScore = A_update_Modularity;

	return(mat);


}

/*
 * void Compute_AOfGroup(spmat* A,spmat* Ag,double* s,int sizeofgroup,int flag)
 *
 *compute sub matrix of A into Ag according to s with size sizeofgroup,
 *compute s.t A[i][j] exits in Ag if s[i]==flag and s[j]==flag.
 */
void Compute_AOfGroup(spmat* A,spmat* Ag,double* s,int sizeofgroup,int flag){
	double *sPtr;
	NodeOfRow** lst = ((ListOfNode*)A->private)->lst;
	int* entry = NULL;
	int j, i, n;
	int lengthRowOfSubMatrix, num_col, colOfSubMatrix, sizeofrow;
	int* row;
	int* newrow;
	int *RowPtr, *NewRowPtr;

	n = A->n;
	row = malloc(sizeof(int)*sizeofgroup);
	if(sizeofgroup != 0 && row ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	sPtr = s;
	for(j = 0; j < sizeofgroup; j++){
		while(*sPtr != flag){
			sPtr++;
			lst++;
		}
		entry = (*lst)->row;
		sizeofrow = (*lst)->num_of_nnz;
		RowPtr = row;
		colOfSubMatrix = 0;
		i = 0;
		lengthRowOfSubMatrix = 0;

		while(sizeofrow != 0 && colOfSubMatrix < sizeofgroup && i < n){
			num_col = *entry;
			if(s[i] == flag && num_col == i){
				entry++;
				*RowPtr = colOfSubMatrix;
				colOfSubMatrix++;
				lengthRowOfSubMatrix++;
				RowPtr++;
				sizeofrow--;
				if (0 != sizeofrow)
					num_col = *entry;
			}
			else if(s[i] == flag){
				colOfSubMatrix++;
			}
			if(0 != sizeofrow && num_col <= i){
				entry++;
				sizeofrow--;
			}
			i++;

		}

		newrow = malloc(sizeof(int)*lengthRowOfSubMatrix);
		if(lengthRowOfSubMatrix != 0 && newrow ==NULL){
			printf("A fatal memory error has been detected.. exiting.\n");
			exit(-1);
		}
		RowPtr = row;
		NewRowPtr = newrow;

		for(i = 0; i < lengthRowOfSubMatrix; i++){
			*NewRowPtr = *RowPtr;
			RowPtr++;
			NewRowPtr++;
		}


		Ag->add_row(Ag, newrow, j, lengthRowOfSubMatrix);
		lst++;
		sPtr++;
	}

	free(row);

}

/*
 * spmat* allocate_sub_spmat(spmat* A,double* s,int n,int flag)
 *
 * allocates sub matrix of A of size n,
 * s.t include each A(i,j) s.t s[i]==flag and s[j]==flag.
 */
spmat* allocate_sub_spmat(spmat* A,double* s,int n,int flag){
	spmat* mat = malloc(sizeof(spmat));
	ListOfNode* list = malloc(sizeof(ListOfNode));

	if(mat == NULL || list ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	list->lst = malloc(sizeof(NodeOfRow) * n);
	if(n != 0 && list->lst ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	mat->n = n;
	mat->private = list;
	mat->add_row = add_row;
	mat->free = mat_free;
	mat->mult = mult;
	mat->shift_amount = 0;
	mat->shift_mat = shift_mat;
	mat->Sum_OF_Rows = Sum_OF_Rows;
	mat->AupdateScore = A_update_Modularity;

	Compute_AOfGroup(A, mat, s, n, flag);
	return mat;
}


