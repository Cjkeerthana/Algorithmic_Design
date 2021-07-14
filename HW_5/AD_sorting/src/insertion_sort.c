#include "insertion_sort.h"
#include "total_order.h"
#include "swap.h"

#define ADDR(A, pos, elem_size)    (A + pos*elem_size)

void insertion_sort(void *A, const unsigned int n, 
                    const size_t elem_size, 
                    total_order leq)
{
    for (size_t i = 0; i < n; i++)
    {
        size_t j = i;
        while (j > 0 && leq(ADDR(A, j, elem_size), ADDR(A, (j-1), elem_size)))
        {
            swap(ADDR(A, j, elem_size), ADDR(A, (j-1), elem_size), elem_size);
            j--; 
        }
    }
    
}