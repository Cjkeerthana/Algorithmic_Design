#ifndef __NONSQ_STRASSEN__
#include <stdlib.h>

void nonsq_strassen_matrix_multiplication(float **C, float const *const *const A,
                                            float const *const *const B, size_t m, size_t k, size_t n);
#endif //__NONSQ_STRASSEN__