#ifndef __TEST_NONSQ__
#include "stdlib.h"

double test_nonsq(void (*f)(float **,
	                  float const *const *const,
	                  float const *const *const,
	                  size_t, size_t, size_t), 
	                  float **C, float** A, float **B, size_t m, size_t k, size_t n);

#endif // __TEST_NONSQ__