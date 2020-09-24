#include <binheap.h>
#include <string.h>
#include <stdio.h>

#define PARENT(node)    ((node-1)/2)
#define LEFT_CHILD(node)    (2*node + 1)
#define RIGHT_CHILD(node)   (2*(node+1))

#define VALID_NODE(H, node) ((H)->num_of_elem>(node))

#define ADDR(H, key_pos)   ((H)->A + (key_pos)*(H)->key_size)
#define INDEX_OF(H, addr)   (((addr) - ((H)->A))/((H)->key_size))

#define NODE_POS(H, node)   (&(H)->node_position[node])
#define KEY_POS(H, node)    (&(H)->key_position[node])

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
    unsigned int kpos = *KEY_POS(H, 0);
    return ADDR(H, kpos);
}

void update_node_pos(binheap_type *H){
    for(unsigned int i = 0; i < H->num_of_elem; i++){
        unsigned int j = H->key_position[i];
        H->node_position[j] = i;
    }
}

void swap_keys(binheap_type *H, unsigned int n_a, unsigned int n_b)
{
    unsigned int *a_pos = KEY_POS(H, n_a);
    unsigned int *b_pos = KEY_POS(H, n_b);
    unsigned int *temp = malloc(sizeof(unsigned int));

    memcpy(temp, a_pos,   sizeof(unsigned int));
    memcpy(a_pos, b_pos,   sizeof(unsigned int));
    memcpy(b_pos, temp, sizeof(unsigned int));

    update_node_pos(H);
    free(temp);
}

void heapify(binheap_type *H, unsigned int node)
{
    unsigned int dst_node=node, child;

    do{
        node = dst_node;

        unsigned int parent_kpos = *KEY_POS(H, node);

        //find the minimum among the node & its children
        child = RIGHT_CHILD(node);
        if(VALID_NODE(H, child)){
                unsigned int child_kpos = *KEY_POS(H, child);

                if(H->leq(ADDR(H, child_kpos), ADDR(H, parent_kpos))){
                    dst_node = child;
                }
        }
        child = LEFT_CHILD(node);
        if(VALID_NODE(H, child)){
                unsigned int child_kpos = *KEY_POS(H, child);

                if(H->leq(ADDR(H, child_kpos), ADDR(H, parent_kpos))){
                    dst_node = child;
        }
        // if the minimum is not in node swap the keys
        if (dst_node !=  node){
            swap_keys(H, dst_node, node);
        }
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
                        total_order_type leq)
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
                         const unsigned int max_size,  
                         const size_t key_size, 
                         total_order_type leq)
{
    binheap_type *H = (binheap_type *)malloc(sizeof(binheap_type));

    H->A = A;
    H->num_of_elem = num_of_elem;
    H->max_size = max_size;
    H->key_size = key_size;
    H->leq = leq;
    H->max_order_value = malloc(key_size);
    H->key_position = (unsigned int *)malloc(sizeof(unsigned int)*max_size);
    H->node_position = (unsigned int *)malloc(sizeof(unsigned int)*max_size);
    H->key_pos_size = sizeof(unsigned int);
    if(num_of_elem == 0){
        return H;
    }

    for (unsigned int i = 0; i < num_of_elem; i++){
        H->key_position[i] = i;
        H->node_position[i] = i;
    }

    const void *value = find_the_max(A, num_of_elem, 
                                        key_size, leq);
    memcpy(H->max_order_value, value, key_size);

    //fix the heap property from the second last level
    for(unsigned int i = num_of_elem/2; i>0; i--){
        heapify(H, i);
    }
    heapify(H, 0);
    update_node_pos(H);
    return H;
}

void delete_heap(binheap_type *H)
{
    free(H->max_order_value);
    free(H);
}

const void *decrease_key(binheap_type *H, void *node, const void *value)
{
    
    unsigned int node_idx = INDEX_OF(H, node);
    unsigned int node_kpos = *KEY_POS(H, node_idx);

    //if the node is not valid or the *value > *node, return NULL
    if(!VALID_NODE(H, node_idx) || !(H->leq(value, node))){
        return NULL;
    }

    memcpy(node, value, H->key_size);

    if(H->num_of_elem > 1){
        unsigned int parent_idx = PARENT(node_idx);
        unsigned int parent_kpos = *KEY_POS(H, parent_idx);
        void *parent = ADDR(H, parent_kpos);

        while((node_idx!=0) && (!H->leq(parent, node))){
            swap_keys(H, parent_kpos, node_kpos);

            node = parent;
            node_idx = parent_idx;
            node_kpos = parent_kpos;

            if(parent_idx != 0){
                parent_idx = PARENT(node_idx);
                parent_kpos = *KEY_POS(H, parent_idx);
                parent = ADDR(H, parent_kpos);
            }
        }
    }
    update_node_pos(H);
    return node;
}

const void *insert_value(binheap_type *H, const void *value)
{   
    //if the heap is already full
    if (H->max_size == H->num_of_elem){
        return NULL;
    }
    
    //if the new value > *max_order_value
    if(H->num_of_elem == 0 || !(H->leq(value, H->max_order_value))){
        memcpy(H->max_order_value, value, H->key_size);
    }

    //get the position of the new node
    void *new_node_addr = ADDR(H, H->num_of_elem);
    unsigned int new_node_pos = INDEX_OF(H, new_node_addr);
    memcpy(new_node_addr, H->max_order_value, H->key_size);
    memcpy(KEY_POS(H,H->num_of_elem), &new_node_pos, sizeof(unsigned int));

    //increase the size of the heap
    H->num_of_elem++;

    //decrease the key of the new node
    return decrease_key(H, new_node_addr, value);
}

void print_heap(const binheap_type *H, 
                void (*key_printer)(const void *value))
{
    unsigned int next_level_node = 1; // store the index of the 
                                      // left-most node of the 
                                      // next level
    for(unsigned int node = 0; node < H->num_of_elem; node++){
        if(node == next_level_node){
            printf("\n");
            next_level_node = LEFT_CHILD(node);
        }else{
            printf("\t");
        }
        unsigned int kpos = *KEY_POS(H,node);
        key_printer(ADDR(H, kpos));
    }
    printf("\n");
}