#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include"SetOfGroups.h"
#include "modularity.h"
#include "BHatMatrix.h"
#include "SPBufferset.h"


/**
 *
 *  Project done By: 208145268 and 209105790.
 *  with love..
 *
 *
 * int main(int argc, char **argv)
 *
 * input: input filename, and output filename.
 * output: returns 0 when successfully ran, -1 otherwise.
 */

int main(int argc, char **argv) {
	FILE* input;
	FILE* output;
	BHat* B;

	groups *O = NULL;


	if(argc != 3 || argv[0] == NULL){return -1;}

	srand(time(NULL));
	SP_BUFF_SET();
	input = fopen(argv[1], "rb");
	if(input == NULL){
		printf("Reading from input file has failed.. exiting.\n");
		return -1;
	}

	output = fopen(argv[2], "wb");
	if(output == NULL){
		printf("Writing to output file has failed.. exiting.\n");
		return -1;
	}



	/*Build BHat of graph*/
	B = BHat_compute(input);

	/*Devises the nodes of graph to groups that given a max modularity */
	O = algorithm3(B);

	/*export the groups in Set O to output according to pdf */
	O->exportGroups(O, output);

	O->free(O);

	return 0;
}

