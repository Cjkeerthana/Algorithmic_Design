#include <stdlib.h>
#include "weighted_graph.h"

#define INFINITY 9999999

GRAPH build_graph(NODE* V, int** W, unsigned int size)
{
    GRAPH G;
    G.V = V;
    G.W = W;
    G.size = size;

    return G;
}

void delete_graph(GRAPH* G)
{
    free(G->V);
    free(G->W);
}