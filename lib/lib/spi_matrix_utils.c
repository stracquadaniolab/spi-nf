#include "spi_matrix_utils.h"

void compute_stochastic_matrix(double * pmat, size_t pmat_dim,
                              int * genotype, size_t genotype_len,
                              double * segment_size, size_t segment_len,
                              double nu, double b
                              )
{
  int i, j, k;
  double genotype_size = 0;
  double d_ij = 0.0, d_ji = 0.0, prob = 0.0, prob_sum = 0.0, bsq = b*b;

  // computing genotype len in bp
  for(i=0; i < genotype_len; i++){
    genotype_size += segment_size[abs(genotype[i])-1];
  }

  //
  for (i=0; i < genotype_len; i++){
    for(j=i+1; j < genotype_len; j++){

      // computing the distance first
      d_ij = 0.0;
      for(k=i; k < j; k++){
        d_ij = d_ij + segment_size[abs(genotype[k])-1];
      }
      // computing reciprocate distance
      d_ji = genotype_size - d_ij;

      // computing probability using self avoiding walk formula
      prob = pow(d_ij * d_ji, -3.0 * nu);
      prob = prob * exp(-2.0 * bsq * (pow(d_ij,-2) + pow(d_ji,-2)));

      // setting the probability for i->j and j->i
      pmat[i*pmat_dim + j] = prob;
      pmat[j*pmat_dim + i] = prob;

      // keeping the total sum of the prob
      prob_sum += prob;
    }
  }

  // normalizing probabilities
  prob_sum = 2.0 * prob_sum;
  for (i=0; i < pmat_dim; i++){
    for(j=i+1; j < pmat_dim; j++){
        pmat[i*pmat_dim + j] = exp(log(pmat[i*pmat_dim + j]) - log(prob_sum));
        pmat[j*pmat_dim + i] = exp(log(pmat[j*pmat_dim + i]) - log(prob_sum));
    }
  }
}
