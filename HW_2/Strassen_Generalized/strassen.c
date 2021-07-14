#include "matrix.h"
#include "stdio.h"

void naive_aux(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t m, const size_t k, const size_t n)
{
    for (size_t y = 0; y < m; y++){
    for (size_t x = 0; x < n; x++){
      float value = 0.0;
      for (size_t z = 0; z < k; z++){
        value += A[y + A_f_row][z + A_f_col]*B[z + B_f_row][x + B_f_col];
      }
      C[y + C_f_row][x + C_f_col] = value;
    }
  }
}                  

void sub_matrix_blocks(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t m, const size_t n)
{
    for (size_t y = 0; y < m; y++){
        for(size_t x = 0; x < n; x++){
            C[y + C_f_row][x + C_f_col] = 
                    A[y + A_f_row][x + A_f_col] - B[y + B_f_row][x + B_f_col];
        }
    }
}                    

void sum_matrix_blocks(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t m, const size_t n)
{
    for (size_t y = 0; y < m; y++){
        for(size_t x = 0; x < n; x++){
            C[y + C_f_row][x + C_f_col] = 
                    A[y + A_f_row][x + A_f_col] + B[y + B_f_row][x + B_f_col];
        }
    }
} 

void fix_up(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t n, const size_t m)
{   
    const size_t n1 = n-m;
    // C11 = C + a12*b21 
    float **c11 = allocate_matrix(n1, n1);
    
    naive_aux(c11, A, B,
                 0, 0,
                A_f_row, A_f_col + n1,
                B_f_row + n1, B_f_col,
                n1, m, n1);  
    
    sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)c11,
                        C_f_row, C_f_col,
                        C_f_row, C_f_col,
                        0, 0,
                        n1, n1);
    deallocate_matrix(c11, n1);
    
    // C12 = A11*b12 + a12*b22 
    float **c12 =  allocate_matrix(n1, m);
    float **c121 = allocate_matrix(n1, m);
    naive_aux(c12, A, B,
                    0, 0,
                    A_f_row, A_f_col,
                    B_f_row, B_f_col + n1,
                    n1, n1, m);   
    naive_aux(c121, A, B,
                    0, 0,
                    A_f_row, A_f_col + n1,
                    B_f_row + n-1, B_f_col + n1,
                    n1, m, m);  

    sum_matrix_blocks(C, (const float* const *const)c12, (const float* const *const)c121,
                        C_f_row, C_f_col + n1,
                        0, 0,
                        0, 0,
                        n1, m);     

    deallocate_matrix(c12, n1);
    deallocate_matrix(c121, n1); 
    
    //C21 = a21*B11 + a22*b21 
    float **c21 =  allocate_matrix(m, n1);
    float **c211 = allocate_matrix(m, n1);
    naive_aux(c21, A, B,
                    0, 0,
                    A_f_row + n1, A_f_col,
                    B_f_row, B_f_col,
                    m, n1, n1);   
                    
    naive_aux(c211, A, B,
                    0, 0,
                    A_f_row + n1, A_f_col + n1,
                    B_f_row + n1, B_f_col,
                    m, m, n1);  
    
    sum_matrix_blocks(C, (const float* const *const)c21, (const float* const *const)c211,
                        C_f_row + n1, C_f_col,
                        0, 0,
                        0, 0,
                        m, n1);     
    
    deallocate_matrix(c21, m);
    deallocate_matrix(c211, m);         
    
    //C22 = a21*b12 + a22*b22
    float **c22 =  allocate_matrix(m, m);
    float **c221 = allocate_matrix(m, m);
    naive_aux(c22, A, B,
                    0, 0,
                    A_f_row + n-1, A_f_col,
                    B_f_row, B_f_col + n1,
                    m, n1, m);
    naive_aux(c221, A, B,
                    0, 0,
                    A_f_row + n1, A_f_col + n1,
                    B_f_row + n1, B_f_col + n1,
                    m, m, m);   
    sum_matrix_blocks(C, (const float* const *const)c22, (const float* const *const)c221,
                        C_f_row + n1, C_f_col + n1,
                        0, 0,
                        0, 0,
                        m, m);    
    deallocate_matrix(c22, m);
    deallocate_matrix(c221, m);  
                                 
}
                
void strassen_aux(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t n)
{
    if(n <= 128){
        naive_aux(C,A,B,
                    C_f_row,C_f_col,
                    A_f_row,A_f_col,
                    B_f_row,B_f_col,
                     n, n, n);
        return;
    }
    //printf("inside strassen aux \n");
    //size of the blocks
    size_t n2;
    if (n%2 != 0){
        n2 = (n-1)/2;
    }
    else{
        n2 = n/2; 
    }

    /*float ***S = (float ***)malloc(sizeof(float **) * 10);
    for (size_t i = 0; i < 10; i++){
        S[i] = allocate_matrix(n2, n2);
    }*/

    float ***P = (float ***)malloc(sizeof(float **) * 7);
    for (size_t i = 0; i < 7; i++){
        P[i] = allocate_matrix(n2, n2);
    }

    // S1 = B12 - B22
    float **S1 = allocate_matrix(n2, n2);
    sub_matrix_blocks(S1, B, B, 
                            0, 0,
                            B_f_row, B_f_col + n2,
                            B_f_row + n2, B_f_col + n2,
                            n2, n2);
    
    //P1 = A11 * S1
    //float **P1 = allocate_matrix(n2, n2);
    strassen_aux(P[0], A, (const float* const *const)S1,
                    0, 0,
                    A_f_row, A_f_col,
                    0, 0,
                    n2);
    deallocate_matrix(S1, n2);

    //S2 = A11 + A12
    float **S2 = allocate_matrix(n2, n2);  
    sum_matrix_blocks(S2, A, A, 
                            0, 0,
                            A_f_row, A_f_col,
                            A_f_row, A_f_col + n2,
                            n2, n2);

    //P2 = S2 * B22
    strassen_aux(P[1], (const float* const *const)S2, B,
                    0, 0,
                    0, 0,
                    B_f_row + n2, B_f_col + n2,
                    n2);
    deallocate_matrix(S2, n2);

    //S3 = A21 + A22
    float **S3 = allocate_matrix(n2, n2);
    sum_matrix_blocks(S3, A, A, 
                            0, 0,
                            A_f_row + n2, A_f_col,
                            A_f_row + n2, A_f_col + n2,
                            n2, n2);

    // P3 = S3 * B11 
    //float **P3 = allocate_matrix(n2, n2); 
    strassen_aux(P[2], (const float* const *const)S3, B,
                    0, 0,
                    0, 0,
                    B_f_row, B_f_col,
                    n2);  
    deallocate_matrix(S3, n2);

    // S4 = B21 - B11
    float **S4 = allocate_matrix(n2, n2);
    sub_matrix_blocks(S4, B, B, 
                            0, 0,
                            B_f_row + n2, B_f_col,
                            B_f_row, B_f_col,
                            n2, n2);
    // P4 = A22 * S4
    strassen_aux(P[3], A, (const float* const *const)S4,
                    0, 0,
                    A_f_row + n2, A_f_col + n2,
                    0, 0,
                    n2);
    deallocate_matrix(S4, n2);

    // S5 = A11 + A22 
    float **S5 = allocate_matrix(n2, n2);  
    sum_matrix_blocks(S5, A, A, 
                            0, 0,
                            A_f_row, A_f_col,
                            A_f_row + n2, A_f_col + n2,
                            n2, n2); 

    // S6 = B11 + B22 
    float **S6 = allocate_matrix(n2, n2); 
    sum_matrix_blocks(S6, B, B, 
                            0, 0,
                            B_f_row, B_f_col,
                            B_f_row + n2, B_f_col + n2,
                            n2, n2); 

    // P5 = S5 * S6
     strassen_aux(P[4], (const float* const *const)S5, 
                    (const float* const *const)S6,
                    0, 0,
                    0, 0,
                    0, 0,
                    n2);   
    deallocate_matrix(S5, n2);
    deallocate_matrix(S6, n2);

    //S7 = A12 - A22
    float **S7 = allocate_matrix(n2, n2); 
    sub_matrix_blocks(S7, A, A, 
                            0, 0,
                            A_f_row, A_f_col + n2,
                            A_f_row + n2, A_f_col + n2,
                            n2, n2);      
    // S8 = B21 + B22
    float **S8 = allocate_matrix(n2, n2); 
    sum_matrix_blocks(S8, B, B, 
                            0, 0,
                            B_f_row + n2, B_f_col,
                            B_f_row + n2, B_f_col + n2,
                            n2, n2); 

    //P6 = S7 * S8  
    strassen_aux(P[5], (const float* const *const)S7, 
                    (const float* const *const)S8,
                    0, 0,
                    0, 0,
                    0, 0,
                    n2);    
    deallocate_matrix(S7, n2);
    deallocate_matrix(S8, n2);

    //S9 = A11 - A21
    float **S9 = allocate_matrix(n2, n2); 
    sub_matrix_blocks(S9, A, A, 
                            0, 0,
                            A_f_row, A_f_col,
                            A_f_row + n2, A_f_col,
                            n2, n2);      

    // S10 = B11 + B12
    float **S10 = allocate_matrix(n2, n2); 
    sum_matrix_blocks(S10, B, B, 
                            0, 0,
                            B_f_row, B_f_col,
                            B_f_row, B_f_col + n2,
                            n2, n2); 

    //P7 = S9 * S10  
    //float **P7 = allocate_matrix(n2, n2);
    strassen_aux(P[6], (const float* const *const)S9, 
                    (const float* const *const)S10,
                    0, 0,
                    0, 0,
                    0, 0,
                    n2);   
    deallocate_matrix(S9, n2);
    deallocate_matrix(S10, n2);

    // C11 = P5 + P4 - P2 + P6
    sum_matrix_blocks(C, (const float* const *const)P[4],
                        (const float* const *const)P[3],
                        C_f_row, C_f_col,
                        0, 0,
                        0, 0,
                        n2, n2);
    sum_matrix_blocks(C, (const float* const *const)C,
                        (const float* const *const)P[5],
                        C_f_row, C_f_col,
                        C_f_row, C_f_col,
                        0, 0,
                        n2, n2);
    sub_matrix_blocks(C, (const float* const *const)C,
                        (const float* const *const)P[1],
                        C_f_row, C_f_col,
                        C_f_row, C_f_col,
                        0, 0,
                        n2, n2);  

    // C12 = P1 + P2
    sum_matrix_blocks(C, (const float* const *const)P[0],
                        (const float* const *const)P[1],
                        C_f_row,C_f_col + n2,
                        0, 0,
                        0, 0,
                        n2, n2);   

    // C21 = P3 + P4 
    sum_matrix_blocks(C, (const float* const *const)P[2],
                        (const float* const *const)P[3],
                        C_f_row + n2,C_f_col,
                        0, 0,
                        0, 0,
                        n2, n2);  

    // C22 = P5 + P1 - P3 - P7
    sum_matrix_blocks(C, (const float* const *const)P[4],
                        (const float* const *const)P[0],
                        C_f_row + n2,C_f_col + n2,
                        0, 0,
                        0, 0,
                        n2, n2);  
    sub_matrix_blocks(C, (const float* const *const)C,
                        (const float* const *const)P[2],
                        C_f_row + n2,C_f_col + n2,
                        C_f_row + n2,C_f_col + n2,
                        0,0,
                        n2, n2);
    sub_matrix_blocks(C, (const float* const *const)C,
                        (const float* const *const)P[6],
                        C_f_row + n2,C_f_col + n2,
                        C_f_row + n2,C_f_col + n2,
                        0,0,
                        n2, n2); 

    /*for(size_t i = 0; i < 10; i++){
        deallocate_matrix(S[i], n2);
    }  
    free(S);*/

    for(size_t i = 0; i < 7; i++){
        deallocate_matrix(P[i], n2);
    }
    free(P); 
    /*deallocate_matrix(P1, n2);
    deallocate_matrix(P3, n2);
    deallocate_matrix(P7, n2); */ 

    if(n%2 != 0){
        fix_up(C,  A, B,
                    C_f_row, C_f_col,
                    A_f_row, A_f_col,
                    B_f_row, B_f_col,
                    n,1);
    }                                                                                                                                             
}
void strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n) 
{
 
    strassen_aux(C, A, B,
                0, 0,
                0, 0,
                0, 0, 
                n);
}

