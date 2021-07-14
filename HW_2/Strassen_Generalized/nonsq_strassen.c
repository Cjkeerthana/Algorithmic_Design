#include "matrix.h"
#include "strassen.h"
#include "stdio.h"

void nonsq_strassen_aux(float **C, float const *const *const A,
                            float const *const *const B,
                            const size_t C_f_row, const size_t C_f_col,
                            const size_t A_f_row, const size_t A_f_col,
                            const size_t B_f_row, const size_t B_f_col,
                            const size_t m, const size_t k, const size_t n)
{
    //printf("inside nonsq_strassen\n");
    if (m == k && k == n){
        strassen_aux(C, A, B,
                        C_f_row, C_f_col, 
                        A_f_row, A_f_col,
                        B_f_row, B_f_col,
                        n);
        return; 
    }
    
    size_t minDim = (m < k)?((m < n)? m : n):((k < n)? k : n);               // choose the minimum of three dimensions to make partitions

    if (minDim <= 128){
        naive_aux(C, A, B,
                    C_f_row, C_f_col,
                    A_f_row, A_f_col, 
                    B_f_row, B_f_col,
                    m, k, n);
        return;
    }
    size_t m1 = m-minDim;
    size_t k1 = k-minDim;
    size_t n1 = n-minDim;
    
    // C11 = a11 * b11 if minDim == k
    // C11 = a11 * b11 + a12 * b21 if minDim == m or n
    //printf("calculating a11b11\n");
    strassen_aux(C, A, B,
                    C_f_row, C_f_col, 
                    A_f_row, A_f_col,
                    B_f_row, B_f_col,
                    minDim); 
    if ((minDim == m && m != k) || (minDim == n && n != k)){
        float **a12b21 = allocate_matrix(minDim,minDim);
        //printf("calculating a12b21\n");
        nonsq_strassen_aux(a12b21, A, B,
                            0, 0,
                            A_f_row, A_f_col + minDim,
                            B_f_row + minDim, B_f_col,
                            minDim,k1,minDim);
        sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)a12b21,
                            C_f_row, C_f_col,
                            C_f_row, C_f_col,
                            0, 0,
                            minDim, minDim);
        deallocate_matrix(a12b21, minDim);
        if(m == n){
            return;
        }
    }

    if((minDim == m && m != n) || (minDim == k && k != n)){
        //C12 = a11 * b12 if minDim == k or m
        //C12 = a11 * b12 + a12 * b22 if minDim == m & m != k
        //printf("calculating a11b12\n");
       
        nonsq_strassen_aux(C, A, B,
                            C_f_row, C_f_col + minDim, 
                            A_f_row, A_f_col,
                            B_f_row, B_f_col + minDim,
                            minDim, minDim, n1); 
        if(m == k){
            return;
        }
        if (minDim == m){            
            float **a12b22 = allocate_matrix(minDim,n1);
            //printf("calculating a12b22\n");
            nonsq_strassen_aux(a12b22, A, B,
                                0, 0, 
                                A_f_row, A_f_col + minDim,
                                B_f_row + minDim, B_f_col + minDim,
                                minDim, k1, n1); 
            sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)a12b22,
                                    C_f_row, C_f_col + minDim,
                                    C_f_row, C_f_col + minDim,
                                    0, 0,
                                    minDim, n1);
            deallocate_matrix(a12b22,minDim);
        }                                                                  
    }
    if((minDim == n && n != m)|| (minDim == k && k != m)){
        //C21 = a21 * b11 if minDim == k
        //C21 = a21 * b11 + a22 * b21 if minDim == n
        //printf("calculating a21b11\n");
        nonsq_strassen_aux(C, A, B,
                            C_f_row + minDim, C_f_col,
                            A_f_row + minDim, A_f_col,
                            B_f_row, B_f_col,
                            m1, minDim, minDim);
        if(n == k){
            return;
        }
        if (minDim == n)
        {
            float **a22b21 = allocate_matrix(m1,minDim);
            //printf("calculating a22b21\n");
            nonsq_strassen_aux(a22b21, A, B,
                                0, 0, 
                                A_f_row + minDim, A_f_col + minDim,
                                B_f_row + minDim, B_f_col,
                                m1, k1, minDim); 
            sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)a22b21,
                                C_f_row + minDim, C_f_col,
                                C_f_row + minDim, C_f_col,
                                0, 0,
                                m1, minDim);
            deallocate_matrix(a22b21,m1);
        }
        
    }
    if(minDim == k && k != m && k != n){
        //C22 == a21 * b12 if minDim == k
        //printf("calculating a21b12\n");
        nonsq_strassen_aux(C, A, B,
                            C_f_row + minDim, C_f_col + minDim,
                            A_f_row + minDim, A_f_col,
                            B_f_row, B_f_col + minDim,
                            m1, minDim, n1);
    }
}
void nonsq_strassen_matrix_multiplication(float **C, float const *const *const A,
                                                float const *const *const B, size_t m, size_t k, size_t n) 
{
    if (m == k && k == n){
        strassen_matrix_multiplication(C, A, B, n);
        return;
    }
    nonsq_strassen_aux(C, A, B,
                        0, 0,
                        0, 0,
                        0, 0,
                        m, k, n);
    /*printf("strassen\n");
    for(size_t x = 0; x < m; x++){
        for(size_t y = 0; y < n; y++){
        printf("%f\t", C[x][y]);
        }
        printf("\n");
    }*/
}