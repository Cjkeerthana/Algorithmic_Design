#include "selection_sort.h"
#include "total_order.h"
#include "swap.h"

#define ADDR(A, pos, elem_size)    (A + pos*elem_size)

void selection_sort(void *A, const unsigned int n, 
                    const size_t elem_size, 
                    total_order leq)
{
    for (unsigned int i = n; i > 0; i--)
    {
        unsigned int max_j = 0;
        for(unsigned int j = 1; j <= i; j++){
            if(leq(ADDR(A, max_j, elem_size), ADDR(A, j, elem_size)))
            {
                max_j = j;
            }
        }
        swap(ADDR(A, i, elem_size), ADDR(A, max_j, elem_size),elem_size);
    }
}