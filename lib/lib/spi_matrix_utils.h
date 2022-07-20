#ifndef __SPI_MATRIX_UTILS__
#define __SPI_MATRIX_UTILS__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void compute_stochastic_matrix(double * pmat, size_t pmat_dim,
                              int * genotype, size_t genotype_len,
                              double * segment_size, size_t segment_len,
                              double nu, double b
                            );
#endif
