#ifndef BHATMATRIX_H_
#define BHATMATRIX_H_


typedef struct _BHatMatrix {

	/* Matrix size (n*n) */
	int	n;

	/*sum of degrees*/
	int M;

	/*shift amount of BHatMatrix*/
	double shift_amount;

	/*degrees of nodes*/
	int* degrees;

	/*sum of rows F*/
	double *F_sumOfRow;

	/*pointer to A graph.*/
	struct _spmat* A;

	/* Frees all resources used by B */
	void (*free)(struct _BHatMatrix *B);

	/* multiplies Bhat by the vector v and writes the result into the vector result (result pre-allocated).*/
	void (*BHat_mult)(struct _BHatMatrix *B, const double *v, double *result);

	/*multiplies B by the vector v and writes the result to the vector result (result pre-allocated).*/
	void (*B_mult)(struct _BHatMatrix *B, const double *v, double *result);

	/* returns eigenvalue of B,
	 * and compute eigenvector with initial vector v,
	 * into result (result pre-allocated).*/
	double (*get_eigenValue)(struct _BHatMatrix *B, double *v, double *result);

	/* updates score s.t for i != index -- score[i] =-4*Di_index*s[i]*sIndex.
	 * if (i == index) score[i] = -score[i].*/
	void (*BupdateModularity)(struct _BHatMatrix *B, double *result, double *s, int sIndex, int index);
} BHat;


/* allocate BHat according to pdf,
 *	s.t. file input represent the graph of nodes.*/
BHat* BHat_compute(FILE* input);

/* multiplies two vectors and returns the sum of the results.*/
double Vector_mult(double *v1, double *v2, int n);

/*divides B into two sub Bhat matrixes of sizes n1 and n2, according to s */
BHat** Compute_subBHat(BHat* B,double* s,int n1,int n2);


#endif /* BHATMATRIX_H_ */
