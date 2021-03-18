
#ifndef SETOFGROUPS_H_
#define SETOFGROUPS_H_
#include "BHatMatrix.h"

typedef struct _Group{

	/*number of nodes in the group */
	int size;

	 /*pointer to the sub Bhat matrix relevant to the group.*/
	BHat* B;

	/*this vector represents the indexes of nodes in the group in the complete graph*/
	int* indexes;

	/* Frees all resources used by g */
	void (*free)(struct _Group *g);
}group;


typedef struct _SetOfGroups {
	/*number of groups in Set; */
	int size;

	/* Adds group g to Set Of groups Set*/
	void (*Insert)(struct _SetOfGroups *Set, struct _Group *g);

	/* removes the first group from set and returns it. */
	group* (*Remove)(struct _SetOfGroups *Set);

	/* Frees all resources used by Set */
	void (*free)(struct _SetOfGroups *Set);

	/* Exports the groups to the output file */
	void (*exportGroups)(struct _SetOfGroups *Set,  FILE* output);

	/* Private field for inner implementation*/
	void *private;

} groups;


/* Allocates a new Set of groups */
groups* SetOfGroups_allocate();

/* Allocates a new Group with the relevant size, indexes and sub Bhat matrix. */
group* group_allocate(int size, int* indexes, BHat* B);




#endif /* SETOFGROUPS_H_ */
