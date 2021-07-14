#include <stdio.h>
#include <string.h>
#include "stdlib.h"
#include "time.h"

#include "dijkstra.h"
#include "weighted_graph.h"

#define INFINITY 9999999
int main()
{
        
    struct timespec start, end; 
    unsigned int s = 30000;
    int **W = (int **)malloc(sizeof(int *) * s);

    for (size_t i = 0; i < s; i++) {
      W[i] = (int *)malloc(sizeof(int) * s);
    }
    
    srand(10);
    for (unsigned int i = 0; i < s; i++) {
      for (unsigned int j = 0; j < s; j++) {
        W[i][j] = INFINITY;
      }
    }
    for (unsigned int i = 0; i < s/2; i++) {
      for (unsigned int j = 0; j < s/2; j++) {
        unsigned int x = rand() % s;
        unsigned int y = rand() % s;
        W[x][y] = rand() % 100;
      }
    }
    
    
    NODE* V1 = (NODE*) malloc(sizeof(NODE)*s);
    NODE* V2 = (NODE*) malloc(sizeof(NODE)*s);

    GRAPH G1 = build_graph(V1, W, s);
    GRAPH G2 = build_graph(V2, W, s);


    printf("\n\n          Time Performance for each version of Dijkstra's Algorithm             \n\n");
    
    printf("size\tArray\t\tHeap\n");
    for(size_t i=0; i<10; i++)
    {
        G1.size = s/(1<<(9-i));  
        G2.size = s/(1<<(9-i));

        printf("%d", G1.size);
        
        clock_gettime(CLOCK_REALTIME, &start);
        DIJKSTRA_ARRAY(&G1, 0);
        clock_gettime(CLOCK_REALTIME, &end);

        printf("\t%lf", (end.tv_sec-start.tv_sec) +
                        (end.tv_nsec-start.tv_nsec)/1E9 );

        clock_gettime(CLOCK_REALTIME, &start);
        DIJKSTRA_HEAP(&G2, 0);
        clock_gettime(CLOCK_REALTIME, &end);

        printf("\t%lf\n", (end.tv_sec-start.tv_sec) +
                        (end.tv_nsec-start.tv_nsec)/1E9 );
        
    }
    printf("\n");
    
    free(W);
    free(G1.V);
    free(G2.V);
    
    return 0;
}