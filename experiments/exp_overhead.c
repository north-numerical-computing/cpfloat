/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * Test overhead of MATLAB mex interface.
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

double timedifference(struct timespec *start, struct timespec *end) {
  return
    (end->tv_sec - start->tv_sec) +
    (double)(end->tv_nsec - start->tv_nsec) * 1e-9;
}

int main() {

  /* Input parameters. */
  int roundingmodes[] = {CPFLOAT_RND_NA, CPFLOAT_RND_NZ, CPFLOAT_RND_NE,
                         CPFLOAT_RND_TP, CPFLOAT_RND_TN, CPFLOAT_RND_TZ,
                         CPFLOAT_RND_SP, CPFLOAT_RND_SE, CPFLOAT_RND_OD};
  size_t nroundings = 9;
  size_t nsizes = 19;
  size_t *sizes = malloc(nsizes * sizeof(size_t));
  size_t i, j, mult = 10;
  for (i = 1; i<3; i++) {
    mult *= 10;
    for (j = 1; j<10; j++)
      sizes[9*(i-1)+(j-1)] = mult*j;
  }
  sizes[9*(i-1)] = mult*10;
  /* for (i=0; i<19; i++) */
  /*   printf("%5lu\n", sizes[i]); */

  size_t ntests = 10;

  optstruct *fpopts = init_optstruct();

  strcpy(fpopts->format, "h");
  fpopts->precision = 11; // t
  fpopts->emax = 15; // emax
  fpopts->subnormal = CPFLOAT_SUBN_USE;
  fpopts->explim = CPFLOAT_EXPRANGE_TARG;
  fpopts->round = CPFLOAT_RND_NE;
  fpopts->flip = CPFLOAT_NO_SOFTERR;
  fpopts->p = 0.5;

  float fmin = ldexpf(1.,-14);

  struct timespec *start = malloc(sizeof(struct timespec));
  struct timespec *end = malloc(sizeof(struct timespec));
  size_t k, n;
  size_t round;
  double *medtimings = malloc(nsizes * nroundings * sizeof(double));
  double *timing = malloc(ntests * sizeof(double));

  /* Binary32 */
  printf("\n*** BINARY32 ***\n");
  const char outfilefloat [] = "./overhead-clang-single.dat";
  FILE *fidf = fopen(outfilefloat, "w");
  float *Xf = malloc(sizes[nsizes-1] * sizes[nsizes-1] * sizeof(float));
  float *Yf;
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    for (j=0; j<n; j++) // generate normal numbers
      Xf[j] = fmin + rand() / (float)RAND_MAX;
    fprintf(fidf, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (round = 0; round < nroundings; round++) {
      fpopts->round = roundingmodes[round];
      for (k=0; k<ntests; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        Yf = malloc(n * sizeof(float));
        cpfloatf(Yf, Xf, n, fpopts);
        free(Yf);
        clock_gettime(CLOCK_MONOTONIC, end);
        timing[k] = timedifference(start, end);
      }
      qsort(timing, ntests, sizeof(*timing), cmpfun);
      medtimings[round + nroundings * i] = timing[ntests/2];
      /* for(k=0; k<ntests; k++) */
      /* fprintf(f, "%.5e  ", timing[k]); */
      fprintf(fidf, " %10.5e", medtimings[round + nroundings * i]);
      fprintf(stdout, " %10.5e", medtimings[round + nroundings * i]);
    }
    fprintf(fidf, "\n");
    fprintf(stdout, "\n");
  }
  free(Xf);
  fclose(fidf);

  /* Binary64 */
  printf("\n*** BINARY64 ***\n");
  const char outfiledouble [] = "./overhead-clang-double.dat";
  FILE *fidd = fopen(outfiledouble, "w");
  double *Xd = malloc(sizes[nsizes-1] * sizes[nsizes-1] * sizeof(double));
  double *Yd;
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    for (j=0; j<n; j++) // generate normal numbers
      Xd[j] = fmin + rand() / (double)RAND_MAX;
    fprintf(fidd, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (round = 0; round < nroundings; round++) {
      fpopts->round = roundingmodes[round];
      for (k=0; k<ntests; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        Yd = malloc(n * sizeof(double));
        cpfloat(Yd, Xd, n, fpopts);
        free(Yd);
        clock_gettime(CLOCK_MONOTONIC, end);
        timing[k] = timedifference(start, end);
      }
      qsort(timing, ntests, sizeof(*timing), cmpfun);
      medtimings[round + nroundings * i] = timing[ntests/2];
      /* for(k=0; k<ntests; k++) */
      /* fprintf(f, "%.5e  ", timing[k]); */
      fprintf(fidd, " %10.5e", medtimings[round + nroundings * i]);
      fprintf(stdout, " %10.5e", medtimings[round + nroundings * i]);
    }
    fprintf(fidd, "\n");
    fprintf(stdout, "\n");
  }
  free(Xd);
  fclose(fidd);

  free(timing);
  free(medtimings);

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
