#include "matrix.h"
#include "strassen.h"
#include "stdio.h"
/*
The algorithm utilizes the strassen multiplication for rectangular matrices
The main idea is to partition the rectangular matrix into square and non square regions
typically like
||      |          ||
||  sq  | non - sq || 
||      |          ||
We utlitize the strassen algorithm to multiply square matrices and recursivley split the 
non-squre regions into sqaure regions using this algorithm.
These partition is done as per the minimum dimension of the rectangular matrices given for the 
mulitiplication.
For example: if A = 2205 * 1450  B = 1450 * 500; the partition is done as 500 * 500 square region
and the remaining as non square region.
Hence the matrix A is splitted as 4 regions 500 * 500, 1705 * 500, 500 * 950, 1705 * 950.
the matrix B is splitted as 500 * 500,  950 * 500.
now we get A as a 2 * 2 matrix with the 4 regions as blocks and B as a 2 * 1 matrix with 2 blocks.
We can now employ a naive multiplication for 2 * 2 and 2 * 1 outer matrix. 
For the multiplying the blocks of the outer matrix, we use strassen if the blocks are square and 
if the blocks are rectangular, we recursively use the nonsq_strassen to create partitions recursively.
*/
void nonsq_strassen_aux(float **C, float const *const *const A,
                            float const *const *const B,
                            const size_t C_f_row, const size_t C_f_col,
                            const size_t A_f_row, const size_t A_f_col,
                            const size_t B_f_row, const size_t B_f_col,
                            const size_t m, const size_t k, const size_t n)
{
    if (m == k && k == n){
        strassen_aux(C, A, B,
                        C_f_row, C_f_col, 
                        A_f_row, A_f_col,
                        B_f_row, B_f_col,
                        n);
        return; 
    }
    
    size_t minDim = (m < k)?((m < n)? m : n):((k < n)? k : n); // choose the minimum of three dimensions to make partitions              

    if (minDim <= 128){                 // if the minimum dimension is less than 128, we opt to naive matrix multiplication
        naive_aux(C, A, B,
                    C_f_row, C_f_col,
                    A_f_row, A_f_col, 
                    B_f_row, B_f_col,
                    m, k, n);
        return;
    }
    //dimensions of the non-square regions
    size_t m1 = m-minDim;       
    size_t k1 = k-minDim;
    size_t n1 = n-minDim;
    
    //Depending upon which dimension is the minimum dimension, the paritioning becomes different.
    //The program is written for each configuration of matrices
    
    // The first block is always square and hence we utilize strassen here.
    // The other blocks are rectangular, where nonsq_strassen is recursively used.

    // C11 = a11 * b11 if minDim == k
    strassen_aux(C, A, B,
                    C_f_row, C_f_col, 
                    A_f_row, A_f_col,
                    B_f_row, B_f_col,
                    minDim);
    // C11 = a11 * b11 + a12 * b21 if minDim == m or n && m != n && n !=k 
    if ((minDim == m && m != k) || (minDim == n && n != k)){
        float **a12b21 = allocate_matrix(minDim,minDim);
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

    //C12 = a11 * b12 if minDim == k or m and m != n and k != n
    if((minDim == m && m != n) || (minDim == k && k != n)){        
        nonsq_strassen_aux(C, A, B,
                            C_f_row, C_f_col + minDim, 
                            A_f_row, A_f_col,
                            B_f_row, B_f_col + minDim,
                            minDim, minDim, n1); 
        if(m == k){
            return;
        }
        //C12 = a11 * b12 + a12 * b22 if minDim == m & m != k
        if (minDim == m){            
            float **a12b22 = allocate_matrix(minDim,n1);
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

    //C21 = a21 * b11 if minDim == k or n and n != m and k != m
    if((minDim == n && n != m)|| (minDim == k && k != m)){
        nonsq_strassen_aux(C, A, B,
                            C_f_row + minDim, C_f_col,
                            A_f_row + minDim, A_f_col,
                            B_f_row, B_f_col,
                            m1, minDim, minDim);
        if(n == k){
            return;
        }
        //C21 = a21 * b11 + a22 * b21 if minDim == n
        if (minDim == n)
        {
            float **a22b21 = allocate_matrix(m1,minDim);
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

    //C22 == a21 * b12 if minDim == k and k != m and k != n
    if(minDim == k && k != m && k != n){
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
}