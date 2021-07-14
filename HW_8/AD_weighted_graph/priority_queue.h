#ifndef __PRIORITY_QUEUE__
#define __PRIORITY_QUEUE__

#include "weighted_graph.h"
#include "stdlib.h"

typedef struct{
    NODE** Qnode;
    unsigned int size;
}queue;

int is_empty(queue* Q);

queue build_queue(GRAPH* G);

NODE* extract_min_queue(queue* Q);

void delete_queue(queue* Q);
#endif