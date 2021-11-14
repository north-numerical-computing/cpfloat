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

void chop_mpfr(mpfr_t *Yd, const double* Xd, size_t n,
               size_t precision, int emax, int rounding_mode) {

  mpfr_set_default_prec(precision);
  mpfr_set_emax(emax);

  for (size_t i = 0; i < n; i++)
    mpfr_set_d(Yd[i], Xd[i], rounding_mode);
}

int main() {

  /* Input parameters. */
  size_t i, j, k, l, n;
  size_t ntests = 10;

  size_t nsizes = 28;
  size_t *sizes = malloc(nsizes * sizeof(double));
  size_t mult = 1;
  for (i = 1; i <= 3; i++) {
    mult *= 10;
    for (j = 1; j < 10; j++)
      sizes[9*(i-1)+(j-1)] = mult*j;
  }
  sizes[9*(i-1)] = mult*10;

  float fmin = ldexpf(1.,-14);

  struct timespec *start = malloc(sizeof(struct timespec));
  struct timespec *end = malloc(sizeof(struct timespec));

  double medtiming;

  /* Binary64 */
  size_t precision [] = {11, 8, 11};
  size_t emax [] = {15, 127, 127};
  size_t nformats = 3;
  printf("\n*** BINARY64 ***\n");
  const char outfile_conv_seq [] = "./timing-mpfr-conv-seq.dat";
  FILE *fids_conv = fopen(outfile_conv_seq, "w");
  const char outfile_op_seq[] = "./timing-mpfr-op-seq.dat";
  FILE *fids_op = fopen(outfile_op_seq, "w");
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    double *X = malloc(n * sizeof(double));
    double *Y = malloc(n * sizeof(double));
    for (j = 0; j < n; j++) { // generate normal numbers
      X[j] = fmin + rand() / (double)RAND_MAX;
      Y[j] = fmin + rand() / (double)RAND_MAX;
    }
    double *timing = malloc(ntests * sizeof(double));
    fprintf(fids_conv, "%6lu", sizes[i]);
    fprintf(fids_op, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (l = 0; l < nformats; l++) {
      /* fpopts->precision = precision[l]; */
      /* fpopts->emax = emax[l]; */

      // Test conversion with allocation.
      for (k = 0; k < ntests; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        mpfr_t *Xd = malloc(n * sizeof(*Xd));
        for (size_t i = 0; i < n; i++)
          mpfr_init(Xd[i]);
        chop_mpfr(Xd, X, n, precision[l], emax[l], MPFR_RNDN);
        clock_gettime(CLOCK_MONOTONIC, end);
        for (size_t i = 0; i < n; i++)
          mpfr_clear(Xd[i]);
        free(Xd);
        timing[k] = timedifference(start, end);
      }
      qsort(timing, ntests, sizeof(*timing), cmpfun);
      medtiming = timing[ntests / 2];
      fprintf(fids_conv, " %10.5e", medtiming);
      fprintf(stdout, " [C] %10.5e", medtiming);

      // Test arithmetic operations.
      mpfr_t *Xd = malloc(n * sizeof(*Xd));
      mpfr_t *Yd = malloc(n * sizeof(*Yd));
      mpfr_t *Zd = malloc(n * sizeof(*Zd));
      for (size_t i = 0; i < n; i++) {
        mpfr_init(Xd[i]);
        mpfr_init(Yd[i]);
        mpfr_init(Zd[i]);
      }
      chop_mpfr(Yd, Y, n, precision[l], emax[l], MPFR_RNDN);
      for (k = 0; k < ntests; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        for (j = 0; j < n; j++)
          mpfr_add(Zd[i], Xd[i], Yd[i], MPFR_RNDN);
        clock_gettime(CLOCK_MONOTONIC, end);
        timing[k] = timedifference(start, end);
      }
      for (size_t i = 0; i < n; i++) {
        mpfr_clear(Xd[i]);
        mpfr_clear(Yd[i]);
        mpfr_clear(Zd[i]);
      }
      free(Xd);
      free(Yd);
      free(Zd);
      qsort(timing, ntests, sizeof(*timing), cmpfun);
      medtiming = timing[ntests/2];
      fprintf(fids_op, " %10.5e", medtiming);
      fprintf(stdout, " [A] %10.5e", medtiming);
      fprintf(stdout, " |");
    }

    free(X);
    free(Y);
    fprintf(fids_conv, "\n");
    fprintf(fids_op, "\n");
    fprintf(stdout, "\n");
  }
  fclose(fids_conv);

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
