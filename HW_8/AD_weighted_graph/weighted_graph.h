#ifndef __WEIGHTED_GRAPH__
#define __WEIGHTED_GRAPH__

#include "stdlib.h"

typedef struct{
    unsigned int no;
    int pred;
    int dist;
}NODE;

typedef struct{
    NODE* V;
    int **W;
    unsigned int size;
}GRAPH;

GRAPH build_graph(NODE* V, int** W, unsigned int size);

void delete_graph(GRAPH* G);
#endif