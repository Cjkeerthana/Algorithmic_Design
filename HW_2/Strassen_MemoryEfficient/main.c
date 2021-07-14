#include <stdio.h>

#include "test.h"
#include "test_nonsq.h"
#include "matrix.h"
#include "strassen.h"
#include  "nonsq_strassen.h"

int main(int argc, char *argv[]) {
  size_t n = 1 << 12;
  size_t m = 5000;
  size_t k = 5000;

  float **A = allocate_random_matrix(n, n);
  float **B = allocate_random_matrix(n, n);
  float **C0 = allocate_matrix(n, n);
  float **C1 = allocate_matrix(n, n);

  printf ("Square matrix multiplication \n");
  printf("n\tStrassen's Alg.\tNaive Alg.\tSame result\n");
  for (size_t j = 1; j <= n; j *= 2) {

    printf("%ld\t", j);
    fflush(stdout);

    printf("%lf\t", test(strassen_matrix_multiplication, C1, A, B, j));
    fflush(stdout);
    printf("%lf\t", test(naive_matrix_multiplication, C0, A, B, j));
    fflush(stdout);
  
    printf("%d\n", same_matrix((float const *const *const)C0,
                               (float const *const *const)C1, j, j));
  }

  deallocate_matrix(A, n);
  deallocate_matrix(B, n);
  deallocate_matrix(C0, n);
  deallocate_matrix(C1, n);
  

  n = 5000;

  float **C = allocate_random_matrix(m, k);
  float **D = allocate_random_matrix(k, n);
  float **C2 = allocate_matrix(m, n);
  float **C3 = allocate_matrix(m, n);

  printf ("Rectangular matrix multiplication \n");

  printf("m\tk\tn\tStrassen's Alg.\tNaive Alg.\tSame result\n");
  for (size_t i = 2344; i <= m; i += 500){
    for (size_t l = 3345; l <= k; l += 500){
      for (size_t j = 1775; j <= n; j += 500) {

        printf("%ld\t%ld\t%ld\t", i, l, j);
        fflush(stdout);

        printf("%lf\t", test_nonsq(nonsq_strassen_matrix_multiplication, C2, C, D, i, l, j));
        fflush(stdout);
        printf("%lf\t", test_nonsq(nonsq_naive_matrix_multiplication, C3, C, D, i, l, j));
        fflush(stdout);
        
        printf("%d\n", same_matrix((float const *const *const)C2,
                                    (float const *const *const)C3, i, j));
      }
    }
  }
 
  deallocate_matrix(C, m);
  deallocate_matrix(D, k);
  deallocate_matrix(C2, m);
  deallocate_matrix(C3, m);

  return 0;
}
