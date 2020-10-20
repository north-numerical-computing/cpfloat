/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * Compare new code with mpfr.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include "timing_fun.h"

void chop_mpfr(double* Yd, const double* Xd, size_t n,
               size_t precision, int emax, int rounding_mode) {

  mpfr_set_default_prec(precision);
  mpfr_set_emax(emax);

  mpfr_t tmp;
  mpfr_init(tmp);

  size_t i;
  for (i=0; i<n; i++) {
    mpfr_set_d(tmp, Xd[i], rounding_mode);
    Yd[i] = mpfr_get_d(tmp, MPFR_RNDN);
  }

}

int main() {

  /* Input parameters. */
  size_t i, j, k, l, n;
  size_t ntests = 10;

  size_t nsizes = 28;
  size_t *sizes = malloc(nsizes*sizeof(double));
  size_t mult = 1;
  for (i = 1; i<=3; i++) {
    mult *= 10;
    for (j = 1; j<10; j++)
      sizes[9*(i-1)+(j-1)] = mult*j;
  }
  sizes[9*(i-1)] = mult*10;

  float fmin = ldexpf(1.,-14);

  struct timespec *start = malloc(sizeof(struct timespec));
  struct timespec *end = malloc(sizeof(struct timespec));

  double *medtimings = malloc(nsizes * sizeof(double));

  /* Binary32 */

  /* Binary64 */
  size_t precision [] = {11, 8, 11};
  size_t emax [] = {15, 127, 127};
  size_t nformats = 3;
  printf("\n*** BINARY64 ***\n");
  const char outfileseq [] = "./timing-mpfr-seq.dat";
  FILE *fids = fopen(outfileseq, "w");
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    double *Xd = malloc(n * sizeof(double));
    double *Yd;
    for (j=0; j<n; j++) // generate normal numbers
      Xd[j] = fmin + rand() / (double)RAND_MAX;
    double *timing = malloc(ntests * sizeof(double));
    fprintf(fids, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (l=0; l<nformats; l++) {
      /* fpopts->precision = precision[l]; */
      /* fpopts->emax = emax[l]; */
      for (k=0; k<ntests; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        Yd = malloc(n * sizeof(double));
        chop_mpfr(Yd, Xd, n, precision[l], emax[l], MPFR_RNDN);
        free(Yd);
        clock_gettime(CLOCK_MONOTONIC, end);
        timing[k] = timedifference(start, end);
      }
      qsort(timing, ntests, sizeof(*timing), cmpfun);
      medtimings[i] = timing[ntests/2];
      fprintf(fids, " %10.5e", medtimings[i]);
      fprintf(stdout, " %10.5e", medtimings[i]);
    }
    free(Xd);
    fprintf(fids, "\n");
    fprintf(stdout, "\n");
  }
  fclose(fids);

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
