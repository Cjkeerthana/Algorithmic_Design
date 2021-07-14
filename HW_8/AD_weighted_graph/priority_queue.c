#include "weighted_graph.h"
#include "priority_queue.h"

int is_empty(queue* Q)
{
    return (Q->size == 0);
}

queue build_queue(GRAPH* G)
{
    queue Q;
    Q.Qnode = (NODE**) malloc(sizeof(NODE*)*G->size);
    for(unsigned int i = 0; i < G->size; i++){
        Q.Qnode[i] = &(G->V[i]);
    }
    Q.size = G->size; 
    return Q;
}

NODE* extract_min_queue(queue* Q)
{
    unsigned int min_index = 0;
    NODE* m = (Q->Qnode[min_index]);
    int min_dist = m->dist;

    for(unsigned int i = 0; i < Q->size; i++){
        NODE* j = (Q->Qnode[i]);
        if(j->dist < min_dist){
            min_dist = j->dist;
            min_index = j->no;
        }
    }

    m = (Q->Qnode[min_index]);
    unsigned int size = Q->size;
    Q->Qnode[min_index] = Q->Qnode[size-1];
    Q->size--;
    return m; 
}

void delete_queue(queue* Q)
{
    free(Q->Qnode);
}
