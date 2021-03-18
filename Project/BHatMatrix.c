
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPBufferset.h"
#include "BHatMatrix.h"
#include "spmat.h"
#define EPSILON 0.00001

/*
 * void B_free(struct _BHatMatrix *B)
 *
 *  Frees all resources used by B
 */
void B_free(struct _BHatMatrix *B){
	free(B->degrees);
	free(B->F_sumOfRow);
	B->A->free(B->A);
	free(B);
}

/*
 * double Vector_mult(double *v1, double *v2, int n)
 *
 *  multiplies two vectors with size n and returns the sum of the results.
 */
double Vector_mult(double *v1, double *v2, int n){
	int i;
	double sum = 0;

	for(i = 0; i < n; i++){
		sum += *v1 * *v2;
		v1++;
		v2++;
	}
	return sum;
}

/*
 * void Degrees_mult(int *degrees,const double *v,double *result, int M, int sizeofgroup)
 *
 * multiplies D by v s.t D[i][j]=-degrees[i]*degrees[j]/M with size sizeofgroup*sizeofgroup,
 * into result (result is pre-allocated).
 */
void Degrees_mult(int *degrees,const double *v,double *result, int M, int sizeofgroup){
	int *K1;
	int i;
	double sum;

	K1 = degrees;
	sum = 0;
	for(i = 0; i < sizeofgroup; i++){
		sum += (double)(*K1) * *v;
		K1++;
		v++;
	}

	K1=degrees;
	for(i = 0; i < sizeofgroup; i++){
		if(*K1!=0)
			*result -= (double)(*K1)*sum/M;
		K1++;
		result++;
	}
}

/*
 * void B_mult(struct _BHatMatrix *B, const double *v, double *result)
 *
 * multiplies matrix B by the vector v,
 * s.t. B[i][i]=A[i][i]+B->F_sumOfRow[i] and B[i][j]=A[i][j] for i!=j,
 * into result (result pre-allocated)
*/
void B_mult(struct _BHatMatrix *A, const double *v, double *result){
	const spmat* mat;

	mat = A->A;
	/*part1--*/
	mat->mult(mat, v, result);

	/*part2--*/
	Degrees_mult(A->degrees, v, result, A->M, A->n);

}

/*
 * void  BHat_mult(struct _BHatMatrix *B, const double *v, double *result)
 *
 * multiplies B by the vector v,
 * into result (result pre-allocated)
*/
void BHat_mult(struct _BHatMatrix *B, const double *v, double *result){
	const spmat* mat;
	int j;
	int n;
	double *FPtr;
	n = B->n;
	mat = B->A;

	/*part1--*/
	mat->mult(mat, v, result);

	/*part2--*/
	Degrees_mult(B->degrees, v, result, B->M, n);


	/*part3--*/
	FPtr = B->F_sumOfRow;

	/*part4--*/
	for(j = 0; j < n; j++){
		*result = *result - (*v * *FPtr);
		FPtr++;
		v++;
		result++;
	}
}

/*
 * double* unitVector(int n)
 *
 * allocates vector with size n s.t vector[i]=1 for 0<=i<=n.
 */
double* unitVector(int n){
	double* v;
	int i;
	double *vPtr;

	v = malloc(sizeof(double)*n);
	if(n != 0 && v ==NULL){
		printf("Error- can't allocate!\n");
		exit(-1);
	}

	vPtr = v;
	for(i = 0; i < n; i++){
		*vPtr = 1;
		vPtr++;
	}
	return v;
}

/*
 * double* F_sumOfRow (const spmat* mat, int M,int* degrees)
 *
 * returns sum of rows of B,
 * s.t. B[i][j]=mat[i][j]-degrees[i]*degrees[j]/M
 */
double* F_sumOfRow (const spmat* mat, int M,int* degrees){
	double* sumOfRows;
	double* unit_vector;
	int n;

	n = mat->n;
	sumOfRows = calloc(n,sizeof(double));
	if(n != 0 && sumOfRows ==NULL){
		printf("Error- can't allocate!\n");
		exit(-1);
	}
	unit_vector = unitVector(n);

	Degrees_mult(degrees, unit_vector, sumOfRows, M, n);
	mat->Sum_OF_Rows(mat, sumOfRows);

	free(unit_vector);
	return sumOfRows;
}

/*
 * void get_eigenVector(struct _BHatMatrix *B, double *v, double *result)
 *
 * compute eigenvector of B by power iteration with initial vector v,
 * into result (result pre-allocated)
 */
void get_eigenVector(struct _BHatMatrix *B, double *v, double *result){
	register int n = 0, i;
	register double dotProduct = 0, distance;
	register double *vPtr, *resultPtr;
	int size;
	int iter = 0;
	int limit;
	size = B->n;

	limit = 200000 + (B->n * 100);
	resultPtr = v;
	for(i = 0; i < size; i++){
		dotProduct += (*resultPtr) * (*resultPtr);
		resultPtr++;
	}
	dotProduct = sqrt(dotProduct);
	resultPtr = v;
	for(i = 0; i < size; i++){
		*resultPtr = (double)(*resultPtr / dotProduct);
		resultPtr++;
	}

	while(n != size){
			iter++;
			n = 0;
			dotProduct = 0;
			B->BHat_mult(B, v, result);

			vPtr = result;
			for(i = 0; i < size; i++){
				dotProduct += (*vPtr) * (*vPtr);
				vPtr++;
			}

			dotProduct = sqrt(dotProduct);
			if(dotProduct == 0){
				printf("Dividsion by zero has been detected.. exiting.\n");
				exit(-1);
			}
			vPtr = result;
			resultPtr = v;
			for(i = 0; i < size; i++){

				*vPtr = (double)(*vPtr / dotProduct);

				distance = (double)fabs(*resultPtr - *vPtr);
				if(distance <= EPSILON)
					n++;
				*resultPtr = *vPtr;
				*vPtr = 0;
				resultPtr++;
				vPtr++;
			}
			if(iter > limit){
				printf("Infinite loop has been detected.. exiting.\n");
				exit(-1);
			}
		}

	vPtr = result;
	resultPtr = v;
	for(i = 0; i < size; i++){
		*vPtr = *resultPtr;
		resultPtr++;
		vPtr++;
	}
}

/*
 * double get_eigenValue(struct _BHatMatrix *B, double *v, double *result )
 *
 * returns eigenvalue of B,
 * and compute eigenvector with initial vector v,
 * into result (result pre-allocated).
 */
double get_eigenValue(struct _BHatMatrix *B, double *v, double *result ){

	double eig = 0;
	spmat* mat;
	int n;

	/*compute eigenVector--*/
	n = B->n;
	mat = B->A;
	mat->shift_amount = B->shift_amount;
	get_eigenVector(B, v, result);

	/*compute eigenValue--*/
	mat->shift_amount = 0;
	B->BHat_mult(B, result, v);
	eig = Vector_mult(v, result, n);
	eig = eig / Vector_mult(result,result,n);

	return eig;
}

/*
 * void B_update_Modularity(struct _BHatMatrix *B, double *result, double *s, int sIndex, int index)
 *
 * updates score s.t score[i] =-4*B[i][index]*s[i]*sIndex for i != index and score[index]=-score[index].
 */
void B_update_Modularity(struct _BHatMatrix *B, double *score, double *s, int sIndex, int index){

	const spmat* mat;
	int M;
	int *degrees;
	register int i;
	int degreeOfIndex;
	double sum;
	int n;

	n = B->n;
	degrees = B->degrees;
	mat = B->A;
	M = B->M;
	degreeOfIndex = degrees[index];

	mat->AupdateScore(mat, score, s, sIndex, index);

	for(i = 0; i < n; i++){
		if( i != index){
			sum = (double)4 * *s * sIndex * *degrees * degreeOfIndex/M;
			*score += sum;
		}
		else{
			*score = -*score;
		}
		degrees++;
		score++;
		s++;
	}
}

/*
 * int* compute_Sub_degrees(int* degrees, double*s, int flag, int n, int n1)
 *
 * returns sub degrees ,s.t:
 * subDegrees with size n1 and degrees with size n,
 * subDegrees include degrees[i] if s[i]==flag.
 */
int* compute_Sub_degrees(int* degrees, double*s, int flag, int n, int n1){
	int* subDegrees;
	int i, j = 0;
	int *subDegreesPtr;

	subDegrees = malloc(sizeof(int) * n1);
	if(n1 != 0 && subDegrees == NULL){
		printf("Error- can't allocate!\n");
	}

	subDegreesPtr=subDegrees;
	for(i = 0; i < n1; i++){
		while(*s != flag && j < n){
			degrees++;
			s++;
			j++;
		}
		*subDegreesPtr = *degrees;
		subDegreesPtr++;
		degrees++;
		s++;
		j++;
	}
	return subDegrees;
}

/*
 * BHat* BHat_allocate(spmat* A, int *degrees, int n, double shift_amount, int M)
 *
 * allocate BHat, s.t:
 * B[i][j]=A[i][j]-degrees[i]*degrees[j]/M for i!=j
 * B[i][j]=A[i][i]-degrees[i]*degrees[i]/M-F_sumOfRow(A,M,degrees)[i].
 */
BHat* BHat_allocate(spmat* A, int *degrees, int n, double shift_amount, int M){

	BHat* B;
	B = malloc(sizeof(BHat));
	if(B == NULL){
		printf("Error- can't allocate!\n");
		exit(-1);
	}

	B->A = A;
	B->F_sumOfRow = F_sumOfRow(A, M, degrees);
	B->M = M;
	B->degrees = degrees;
	B->n = n;
	B->shift_amount = shift_amount;

	B->free = B_free;
	B->get_eigenValue = get_eigenValue;
	B->BupdateModularity = B_update_Modularity;
	B->B_mult = B_mult;
	B->BHat_mult = BHat_mult;

	return B;
}

/*
 * BHat** Compute_subBHat(BHat *B,double* s,int n1,int n2)
 *
 * return array BHats with size tow,
 * which include division of B to Tow subBHat according to division by s,
 * s.t BHats[0] with size n1 and BHats with size n2.
 */
BHat** Compute_subBHat(BHat *B,double* s,int n1,int n2){
	spmat* A1;
	spmat* A2;
	int* degrees1;
	int* degrees2;
	int flag;
	BHat** TOWBHat;

	TOWBHat = malloc (sizeof(BHat*) * 2);
	if(TOWBHat==NULL){
		printf("Error- can't allocate!\n");
		exit(-1);
	}

	flag = 1;
	A1 = allocate_sub_spmat(B->A, s, n1, flag);
	degrees1 = compute_Sub_degrees(B->degrees, s, flag, B->n, n1);

	flag = -1;
	A2 = allocate_sub_spmat(B->A,s,n2,flag);
	degrees2 = compute_Sub_degrees(B->degrees, s, flag, B->n, n2);

	TOWBHat[0] = BHat_allocate(A1, degrees1, n1, B->shift_amount, B->M);
	TOWBHat[1] = BHat_allocate(A2, degrees2, n2, B->shift_amount, B->M);
	return TOWBHat;
}

/*
 *	BHat* BHat_compute(FILE* input)
 *
 *	returns struct BHat according to pdf,
 *	s.t. file input represent the graph of nodes.
 */
BHat* BHat_compute(FILE* input){

	BHat* B = malloc(sizeof(BHat));
	int *degreesPtr, nodes;
	int i, M = 0;
	int* degrees;
	int* indexes;
	spmat *mat;

	if(B ==NULL){
		printf("Error- can't allocate!\n");
		exit(-1);
	}

	fread(&nodes, sizeof(int), 1, input);
	degrees = malloc(sizeof(int) * nodes);
	if(nodes != 0 && degrees ==NULL){
		printf("Error- can't allocate!\n");
		exit(-1);
	}


	/*check if nnz=0*/
	mat = spmat_allocate(nodes);

	degreesPtr = degrees;
	/*
	 * build mat & calculates M.
	*/

	for(i = 0; i < nodes; i++){
		fread(degreesPtr, sizeof(int), 1, input);
		indexes = malloc(sizeof(int) *(*degreesPtr));
		if(*degreesPtr != 0 && indexes ==NULL){
			printf("Error- can't allocate!\n");
			exit(-1);
		}
		fread(indexes, sizeof(int), *degreesPtr, input);


		mat->add_row(mat, indexes, i, *degreesPtr);

		M += *degreesPtr;
		degreesPtr++;
	}

	fclose(input);
	if(M==0){
		printf("Error- can't divide by zero!\n");
		exit(-1);
	}

	B->A = mat;
	B->M = M;
	B->degrees = degrees;
	B->n = nodes;
	B->F_sumOfRow = calloc(nodes,sizeof(double));
	if(nodes != 0 && B->F_sumOfRow ==NULL){
		printf("Error- can't allocate!\n");
		exit(-1);
	}
	B->shift_amount = mat->shift_mat(mat, degrees, M);
	B->free = B_free;
	B->get_eigenValue = get_eigenValue;
	B->BHat_mult = BHat_mult;
	B->BupdateModularity = B_update_Modularity;
	B->B_mult = B_mult;

	return B;
}

