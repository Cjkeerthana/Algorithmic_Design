#include "heap_sort.h"
#include <string.h>
#include <stdlib.h>

#define PARENT(node)    ((node-1)/2)
#define LEFT_CHILD(node)    (2*node + 1)
#define RIGHT_CHILD(node)   (2*(node+1))

#define VALID_NODE(H, node) ((H)->num_of_elem>(node))

#define ADDR(H, node)   ((H)->A+(node)*(H)->key_size)
#define INDEX_OF(H, addr)   (((addr) - ((H)->A))/((H)->key_size))
#define ADDR1(A, pos, elem_size)    (A + pos*elem_size)


typedef struct {
    void *A; // the array used to store heap nodes
    unsigned int num_of_elem; // num of nodes in the heap
    size_t key_size; //size of the key type
    total_order leq; // the total order of the heap
    void *max_order_value; // the maximum value stored in the heap acc. to the order
} binheap_type;

int is_heap_empty(const binheap_type *H)
{
    return H->num_of_elem==0;
}

const void *min_value(const binheap_type *H)
{
    if (is_heap_empty(H)){
        return NULL;
    }
    // the minimum is stored in the root A[0]
    return ADDR(H, 0);
}

void swap_keys(binheap_type *H, unsigned int n_a, unsigned int n_b)
{
    void *p_a = ADDR(H, n_a);
    void *p_b = ADDR(H, n_b);
    void *tmp = malloc(H->key_size);

    memcpy(tmp, p_a, H->key_size);
    memcpy(p_a, p_b, H->key_size);
    memcpy(p_b, tmp, H->key_size);

    free(tmp);
}

void heapify(binheap_type *H, unsigned int node)
{
    unsigned int dst_node=node, child;

    do{
        node = dst_node;

        //find the minimum among the node & its children
        child = RIGHT_CHILD(node);
        if(VALID_NODE(H, child) && 
            H->leq(ADDR(H, child), ADDR(H, dst_node))){
                dst_node = child;
            }
        child = LEFT_CHILD(node);
        if(VALID_NODE(H, child) && 
            H->leq(ADDR(H, child), ADDR(H, dst_node))){
                dst_node = child;
            }
        // if the minimum is not in node swap the keys
        if (dst_node !=  node){
            swap_keys(H, dst_node, node);
        }
    }while(dst_node != node);
}

const void *extract_min(binheap_type *H)
{
    if(is_heap_empty(H)){
        return NULL;
    }

    //swap the keys of root and the  right most leaf of the last level
    swap_keys(H, 0, H->num_of_elem-1);

    //delete the right most leaf of the last level A[num_of_elem-1]
    H->num_of_elem--;
    
    heapify(H, 0);

    return ADDR(H, H->num_of_elem+1);
}

const void *find_the_max(void *A,
                        const unsigned int num_of_elem,
                        const size_t key_size,
                        total_order leq)
{
    if(num_of_elem == 0){
        return NULL;
    }
    const void *max_value = A;
    for(const void *addr = A+key_size; addr!=A+num_of_elem*key_size; 
        addr+=key_size)
    {
        if(leq(addr, max_value)){
            max_value = addr;
        }
    }
    return max_value;
}

binheap_type *build_heap(void *A, 
                         const unsigned int num_of_elem, 
                         const size_t key_size, 
                         total_order leq)
{
    binheap_type *H = (binheap_type *)malloc(sizeof(binheap_type));

    H->A = A;
    H->num_of_elem = num_of_elem;
    H->key_size = key_size;
    H->leq = leq;
    H->max_order_value = malloc(key_size);

    if(num_of_elem == 0){
        return H;
    }

    const void *value = find_the_max(A, num_of_elem, 
                                        key_size, leq);
    memcpy(H->max_order_value, value, key_size);

    //fix the heap property from the second last level
    for(unsigned int i = num_of_elem/2; i>0; i--){
        heapify(H, i);
    }
    heapify(H, 0);

    return H;
}

void delete_heap(binheap_type *H)
{
    free(H->max_order_value);
    free(H);
}

void heap_sort(void *A, const unsigned int n, 
               const size_t elem_size, 
               total_order leq)
{
    binheap_type *H = build_heap(A, n, elem_size, geq_int);

    for (size_t i = n; i > 0; i--)
    {
        const void *b = extract_min(H);
        memcpy(ADDR1(A, i, elem_size), b, elem_size);
    }
    
    delete_heap(H);
}