#ifndef __STRASSEN__

void naive_aux(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t m, const size_t k, const size_t n);

void sum_matrix_blocks(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t m, const size_t n);         

void strassen_aux(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t n);
                    
void strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n);
                                    
#endif //__STRASSEN__
