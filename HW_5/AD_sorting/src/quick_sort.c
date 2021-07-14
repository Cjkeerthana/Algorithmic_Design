#include "quick_sort.h"
#include "swap.h"
#include "total_order.h"

#define ADDR(A, pos, elem_size)    (A + pos*elem_size)

const int partition(void *A, const int l, const int r, const int p, 
                                const size_t elem_size, total_order leq)
{
    int i = l, j = r, k = p;
    swap(ADDR(A, i, elem_size), ADDR(A, k, elem_size), elem_size);
    k = i;
    i = i+1;

    while (i <= j)
    {
        if(leq(ADDR(A, k, elem_size), ADDR(A, i, elem_size))){
            swap(ADDR(A, i, elem_size), ADDR(A, j, elem_size), elem_size);
            j--;
        }
        else
        {
            i++;
        }
    }
    swap(ADDR(A, k, elem_size), ADDR(A, j, elem_size), elem_size);
    return j;
}
void quick_sort_aux(void *A, int l, int r, 
                const size_t elem_size, 
                total_order leq)
{
    while(l < r){
        const int pi = partition(A, l, r, (l+r)/2, elem_size, leq);
        quick_sort_aux(A, l, pi-1, elem_size, leq);
        l = pi + 1;
    }
}
void quick_sort(void *A, const unsigned int n, 
                const size_t elem_size, 
                total_order leq)
{
    quick_sort_aux(A, 0, n-1, elem_size, leq);
}