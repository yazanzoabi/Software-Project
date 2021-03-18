#ifndef SPMAT_H_
#define SPMAT_H_


typedef struct _spmat {
	/* Matrix size (n*n) */
	int		n;

	/* shift amount used to get the maximum eigenvector. */
	double shift_amount;

	/* Adds row i to the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1)
	 */
	void	(*add_row)(struct _spmat *A, int *row, int i,int size);

	/* Frees all resources used by A */
	void	(*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v,
	 * into result (result is pre-allocated).
	 */
	void	(*mult)(const struct _spmat *A, const double *v, double *result);

	/*
	 * calculates the shift amount of a matrix B,
	 * this function is called once, the result is then used for all sub B matrixes.
	 */
	double (*shift_mat)(struct _spmat *A, int *degrees,int M);

	/* ompute sum of each row in A,
	 * into sumOFRows (sumOFRows pre-allocated).
	 */
	void (*Sum_OF_Rows)(const struct _spmat *A,double *sumOFRows);

	/** updates score s.t score[i] =-4*A[i][index]*s[i]*sIndex for i != index.*/
	void (*AupdateScore)(const struct _spmat *A, double *score, double *s, int sIndex, int index);


	/* Private field for inner implementation. */
	void	*private;
} spmat;


/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate(int n);

/* compute sub matrix of A which include each A[i][j] s.t s[i]==flag and s[j]==flag, of size n */
spmat* allocate_sub_spmat(spmat* A,double* s,int n,int flag);

#endif
