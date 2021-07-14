#include <stdlib.h>
#include <stdio.h>
#include "weighted_graph.h"
#include "dijkstra.h"
#include "priority_queue.h"
#include "binheap.h"
#include "total_orders.h"

#define INFINITY 9999999

void INIT_SSSP(GRAPH* G)
{
    for(size_t i = 0; i < G->size; i++){
        NODE* a = &G->V[i];
        a->dist = INFINITY;
        a->pred = INFINITY;
    }
}

void relax(NODE* u, NODE* v, int w)
{
    if(u->dist + w < v->dist){
        v->dist = u->dist + w;
        v->pred = u->no;
    }
}

void DIJKSTRA_ARRAY(GRAPH* G, unsigned int source)
{
    INIT_SSSP(G);
    NODE* src = &G->V[source];
    src->dist = 0;

    queue Q = build_queue(G);
    while(!is_empty(&Q)){
        NODE* u = extract_min_queue(&Q);
        for(unsigned int i = 0; i < G->size; i++){
            int w = G->W[u->no][i];
            if(w < INFINITY){
                NODE* v = &G->V[i];
                relax(u,v,w);
            }
        }
    }
    delete_queue(&Q);
}

void DIJKSTRA_HEAP(GRAPH* G, unsigned int source)
{
    INIT_SSSP(G);
    NODE* src = &G->V[source];
    src->dist = 0;

    binheap_type* H = build_heap(G->V, G->size, G->size, sizeof(NODE), leq_node);

    while(!is_heap_empty(H)){
        NODE* u = (NODE* )extract_min(H);
        for(unsigned int i = 0; i < G->size; i++){
            int w = G->W[u->no][i];
            if(w < INFINITY){
                NODE* v = &G->V[i];
                relax(u,v,w);
                heapify(H,0);
            }
        }
    }
    delete_heap(H);
}