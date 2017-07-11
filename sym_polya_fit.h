#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "seqTM.h"

double sym_polya_fit(double**counts, double* count_sum, int size1, int size2, double initV);
double sym_polya_fit1(int**counts, int* count_sum, int size1, int size2, double initV);
