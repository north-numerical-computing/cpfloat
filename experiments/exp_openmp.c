/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * Compare new code with floatx.
 */

#include <time.h>
#include <stdio.h>
#include <string.h>

#include "cpfloat_binary32.h"
#include "cpfloat_binary64.h"

int cmpfun(const void *x, const void *y) {
  if (*(double *)x < *(double *)y)
    return -1;
  else if (*(double *)x > *(double *)y)
    return 1;
  else
    return 0;
}

double timedifference(struct timespec *start, struct timespec *end)
{
  return
    (end->tv_sec - start->tv_sec) +
    (double)(end->tv_nsec - start->tv_nsec) * 1e-9;
}

int main() {

  /* Input parameters. */
  size_t i, j, k, l, n, r;
  size_t ntests = 100;

  int roundingmodes[] = {CPFLOAT_RND_NA, CPFLOAT_RND_NZ, CPFLOAT_RND_NE,
                         CPFLOAT_RND_TP, CPFLOAT_RND_TN, CPFLOAT_RND_TZ,
                         CPFLOAT_RND_SP, CPFLOAT_RND_SE, CPFLOAT_RND_OD};
  size_t nroundings = 9;

  size_t nsizes = 28;
  size_t *sizes = malloc(nsizes*sizeof(double));
  size_t mult = 1;
  for (i = 1; i <= 3; i++) {
    mult *= 10;
    for (j = 1; j < 10; j++)
      sizes[9*(i-1)+(j-1)] = mult*j;
  }
  sizes[9*(i-1)] = mult*10;

  optstruct *fpopts = init_optstruct();

  fpopts->subnormal = CPFLOAT_SUBN_USE;
  fpopts->explim = CPFLOAT_EXPRANGE_TARG;
  fpopts->flip = CPFLOAT_NO_SOFTERR;
  fpopts->p = 0.5;

  float fmin = ldexpf(1.,-14);

  struct timespec *start = malloc(sizeof(struct timespec));
  struct timespec *end = malloc(sizeof(struct timespec));

  double *medtimings = malloc(nsizes * sizeof(double));

  int numthreads;
  #pragma omp parallel
  {
    numthreads = omp_get_num_threads();
  }
  cpfloat_precision_t precision [] = {11, 8, 11};
  cpfloat_exponent_t emax [] = {15, 127, 127};
  size_t nformats = 3;

  /* Binary32 */

  /* Binary64 */
  printf("\n*** BINARY64 ***\n");
  char * outfiled = malloc(50 * sizeof(outfiled[0]));
  sprintf(outfiled, "./timing-cpfloat-%02d.dat", numthreads);
  char * outfilepard = malloc(50 * sizeof(outfilepard[0]));
  sprintf(outfilepard, "./timing-cpfloat-par-%02d.dat", numthreads);
  const char outfileseqd [] = "./timing-cpfloat-seq.dat";
  FILE *fidd = fopen(outfiled, "w");
  FILE *fidpd = fopen(outfilepard, "w");
  FILE *fidsd = fopen(outfileseqd, "w");
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    double *Xd = malloc(n * sizeof(double));
    double *Yd;
    for (j=0; j<n; j++) // generate normal numbers
      Xd[j] = fmin + rand() / (double)RAND_MAX;
    double *timing = malloc(ntests * sizeof(double));
    fprintf(fidd, "%6lu", sizes[i]);
    fprintf(fidpd, "%6lu", sizes[i]);
    fprintf(fidsd, "%6lu", sizes[i]);
    fprintf(stdout, "*** %6lu\n", sizes[i]);
    for (l=0; l<nformats; l++) {
      fpopts->precision = precision[l];
      fpopts->emax = emax[l];
      for (r = 0; r < nroundings; r++) {
        fpopts->round = roundingmodes[r];

        for (k=0; k<ntests; k++) {
          clock_gettime(CLOCK_MONOTONIC, start);
          Yd = malloc(n * sizeof(double));
          cpfloat(Yd, Xd, n, fpopts);
          free(Yd);
          clock_gettime(CLOCK_MONOTONIC, end);
          timing[k] = timedifference(start, end);
        }
        qsort(timing, ntests, sizeof(*timing), cmpfun);
        medtimings[i] = timing[ntests/2];
        fprintf(fidd, " %10.5e", medtimings[i]);
        fprintf(stdout, " %10.5e", medtimings[i]);

        for (k=0; k<ntests; k++) {
          clock_gettime(CLOCK_MONOTONIC, start);
          Yd = malloc(n * sizeof(double));
          cpfloat_parallel(Yd, Xd, n, fpopts);
          free(Yd);
          clock_gettime(CLOCK_MONOTONIC, end);
          timing[k] = timedifference(start, end);
        }
        qsort(timing, ntests, sizeof(*timing), cmpfun);
        medtimings[i] = timing[ntests/2];
        fprintf(fidpd, " %10.5e", medtimings[i]);
        fprintf(stdout, " %10.5e", medtimings[i]);

        for (k=0; k<ntests; k++) {
          clock_gettime(CLOCK_MONOTONIC, start);
          Yd = malloc(n * sizeof(double));
          cpfloat_sequential(Yd, Xd, n, fpopts);
          free(Yd);
          clock_gettime(CLOCK_MONOTONIC, end);
          timing[k] = timedifference(start, end);
        }
        qsort(timing, ntests, sizeof(*timing), cmpfun);
        medtimings[i] = timing[ntests/2];
        fprintf(fidsd, " %10.5e", medtimings[i]);
        fprintf(stdout, " %10.5e \n", medtimings[i]);
      }
    }
    free(Xd);
    fprintf(fidd, "\n");
    fprintf(fidpd, "\n");
    fprintf(fidsd, "\n");
    fprintf(stdout, "\n");
  }
  fclose(fidd);
  fclose(fidpd);
  fclose(fidsd);

  free_optstruct(fpopts);

  return 0;
}

/*
 * CPFloat - Custom Precision Floating-point numbers.
 *
 * Copyright 2020 Massimiliano Fasi and Mantas Mikaitis
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this library; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */
