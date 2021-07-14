#ifndef __DIJKSTRA__
#define __DIJKSTRA__

#include "stdlib.h"
#include "weighted_graph.h"

GRAPH build_graph(NODE* V, int** W, unsigned int size);

void DIJKSTRA_ARRAY(GRAPH* G, unsigned int source);

void DIJKSTRA_HEAP(GRAPH* G, unsigned int source);

#endif