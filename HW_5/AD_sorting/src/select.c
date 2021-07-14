#include "select.h"
#include "quick_sort.h"
#include "swap.h"
#include "total_order.h"
#include "stdio.h"

#define ADDR(A, pos, elem_size)    (A + pos*elem_size)

int select_alg(void *A, const int l, 
                                const int r, 
                                const int i, 
                                const size_t elem_size, total_order leq);

void partition2(void *A, const int l, const int r, const int p, int* k_1, int* k_2,
                              const size_t elem_size, total_order leq)
{
    int i = l, j = r, k = p, m;
    swap(ADDR(A, i, elem_size), ADDR(A, k, elem_size), elem_size);
    k = i;
    m = k+1;
    i = i+1;

    while (i <= j)
    {
        if(less_int(ADDR(A, k, elem_size), ADDR(A, i, elem_size))){
            swap(ADDR(A, i, elem_size), ADDR(A, j, elem_size), elem_size);
            j--;
        }
        else if(eq_int(ADDR(A, k, elem_size), ADDR(A, i, elem_size))){
            swap(ADDR(A, i, elem_size), ADDR(A, m, elem_size), elem_size);
            m++;
            i++;
        }
        else {
            i++;
        }
    }

    for(int x = k; x < m; x++){
        swap(ADDR(A, x, elem_size), ADDR(A, j, elem_size), elem_size);
        j--;
    }
    *k_1 = j+1;
    *k_2 = j+(m-k); 
}

unsigned int select_pivot(void *A, const int l, const int r, 
                          const size_t elem_size, 
                          total_order leq)
{
    if((r-l+1) <= 10)
    {
        quick_sort_aux(A, l, r, elem_size, leq);
        return (l+r)/2;
    }
    unsigned int chuncks = (r-l+1)/5 ;

    for(int c = 0; c < chuncks; c++)
    {
        unsigned int cl = l + c*5;
        quick_sort_aux(A, cl, cl+4, elem_size, leq);
        swap(ADDR(A, (cl+2), elem_size), ADDR(A, (l+c), elem_size), elem_size);
    }
    
    return select_alg(A, l, (l+chuncks-1), l+(chuncks/2), elem_size, leq);	
}

int select_alg(void *A, const int l, const int r, 
                                const int i, 
                                const size_t elem_size, total_order leq)
{
    if((r-l+1) <= 10)
    {
        quick_sort_aux(A, l, r, elem_size, leq);
        return i;
    }
    unsigned int j = select_pivot(A, l, r, elem_size, leq);
    int k_1, k_2;
    partition2(A, l, r, j, &k_1, &k_2, elem_size, leq);

    if(i < k_1){
        return select_alg(A, l, k_1-1, i, elem_size, leq);
    }
    else if(i > k_2){
        return select_alg(A, k_2+1, r, i, elem_size, leq);
    }
    else{
        return i;
    }
}

void quick_sort_select_aux(void *A, int l, int r, 
                const size_t elem_size, 
                total_order leq)
{    
    while(l < r){
        unsigned int pivot = select_pivot(A, l, r, elem_size, leq);
        int pi_1, pi_2;
        partition2(A, l, r, pivot, &pi_1, &pi_2, elem_size, leq);
        quick_sort_select_aux(A, l, (pi_1-1), elem_size, leq);
        l = pi_2 + 1;
    }
}

void quick_sort_select(void *A, const unsigned int n, 
                       const size_t elem_size, 
                       total_order leq)
{  
   
   quick_sort_select_aux(A, 0, n-1, elem_size, leq);
   
}
