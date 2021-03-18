
#include<stdio.h>
#include<stdlib.h>
#include "SPBufferset.h"
#include "SetOfGroups.h"


/*
 * void group_free(struct _Group *g)
 *
 * Frees all resources used by g
 */
void group_free(struct _Group *g){
	free(g->indexes);
	g->B->free(g->B);
	free(g);
}
/*
 * group* group_allocate(int n,int* indexes,BHat* B)
 *
 * Allocates a new Group with the relevant size, indexes and sub Bhat matrix.
 */
group* group_allocate(int n,int* indexes,BHat* B){
	group* g = malloc(sizeof(group));
	if(g ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}
	g->indexes = indexes;
	g->size = n;
	g->free = group_free;
	g->B=B;
	return g;
}

/*
 * struct _NodeOfGroup
 *
 * represents a node in a group. which includes a pointer to the group and pointer to the next node.
 */
typedef struct _NodeOfGroup {
	group* g;
	struct _NodeOfGroup* next;
} NodeG;

/*
 * struct _ListOfGroups
 *
 * list of of nodesOfGroups and incldes pointers to the first and last nodes in the list.
 */
typedef struct _ListOfGroups {
	NodeG* first;
	NodeG* last;
} ListOfGroups;

/*
 * void groups_insert(struct _SetOfGroups *Set, struct _Group *g)
 *
 * inserts g a new group to a list of groups of Set,
 * and increases the size value by one.
 */
void groups_insert(struct _SetOfGroups *Set, struct _Group *g){
	ListOfGroups* list = NULL;
	NodeG* node = malloc(sizeof(NodeG));

	if(node ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	list = (ListOfGroups*)Set->private;
	node->g = g;
	node->next = NULL;

	if(list->first == NULL){
		list->first = node;
		list->last = node;
	}
	else{
		list->last->next = node;
		list->last = node;
	}

	Set->size++;
}

/*
 * group* groups_remove(struct _SetOfGroups *Set)
 *
 * removes the last added group from Set and returns it.
 */
group* groups_remove(struct _SetOfGroups *Set){
	ListOfGroups* list = NULL;
	group* g = NULL;
	group* newg=NULL;
	NodeG* node = NULL;

	if(Set->size == 0)
		return NULL;
	list = (ListOfGroups*)Set->private;
	g=list->first->g;

	newg=group_allocate(g->size,g->indexes,g->B);
	free(g);
	Set->size--;
	node = list->first;
	list->first = list->first->next;
	free(node);
	return newg;
}

/*
 * void exportGroups(struct _SetOfGroups *Set,  FILE* output)
 *
 * Exports Set into the output file according to instructions from the project's pdf.
 */
void exportGroups(struct _SetOfGroups *Set,  FILE* output){
	int* list;
	int sizeofgroup, size;
	group* g;
	size = Set->size;
	g = Set->Remove(Set);
	if(g == NULL){
		return;
	}
	else{
		fwrite(&size, sizeof(int), 1, output);
		while(g != NULL){
			list = g->indexes;
			sizeofgroup = g->size;
			fwrite(&(g->size), sizeof(int), 1, output);
			fwrite(list, sizeof(int), sizeofgroup, output);
			g->free(g);
			g = Set->Remove(Set);
		}
	}
	fclose(output);
}

/*
 * void groups_free(struct _SetOfGroups *Set)
 *
 * Frees all resources used by Set
 */
void groups_free(struct _SetOfGroups *Set){
	NodeG* node = NULL;
	NodeG* ptr = NULL;
	ListOfGroups* list = NULL;

	list = (ListOfGroups*)Set->private;
	node = list->first;

	while(node != NULL){
		ptr = node;
		ptr->g->free(ptr->g);
		node = node->next;
		free(ptr);
	}
	free(list);
	free(Set);
}

/*
 * groups* SetOfGroups_allocate()
 *
 * Allocates a new Set of groups.
 */
groups* SetOfGroups_allocate(){
	groups* Set = malloc(sizeof(groups));
	ListOfGroups* list = malloc(sizeof(ListOfGroups));

	if(Set == NULL || list ==NULL){
		printf("A fatal memory error has been detected.. exiting.\n");
		exit(-1);
	}

	list->first = NULL;
	list->last = NULL;

	Set->private = list;
	Set->size = 0;
	Set->Insert = groups_insert;
	Set->Remove = groups_remove;
	Set->free = groups_free;
	Set->exportGroups = exportGroups;


	return Set;
}






