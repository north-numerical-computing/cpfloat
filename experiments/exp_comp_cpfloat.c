/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * Compare new code with floatx.
 */

#include <stdio.h>
#include <string.h>

#include "cpfloat_binary32.h"
#include "cpfloat_binary64.h"

#include "timing_fun.h"

int main() {

  /* Input parameters. */
  size_t i, j, k, l, n;

  size_t nsizes = 28;
  size_t *sizes = malloc(nsizes*sizeof(double));
  size_t mult = 1;
  for (i = 1; i <= 3; i++) {
    mult *= 10;
    for (j = 1; j < 10; j++)
      sizes[9 * (i - 1) + (j - 1)] = mult * j;
  }
  sizes[9 * (i - 1)] = mult * 10;

  optstruct *fpopts = init_optstruct();

  strcpy(fpopts->format, "h");
  fpopts->precision = 11; // t
  fpopts->emax = 15;      // emax
  fpopts->subnormal = CPFLOAT_SUBN_USE;
  fpopts->explim = CPFLOAT_EXPRANGE_TARG;
  fpopts->round = CPFLOAT_RND_NE;
  fpopts->flip = CPFLOAT_NO_SOFTERR;
  fpopts->p = 0.5;

  float fmin = ldexpf(1., -14);

  struct timespec *start = malloc(sizeof(struct timespec));
  struct timespec *end = malloc(sizeof(struct timespec));

  double medtiming;

  /* Binary64 */
  cpfloat_precision_t precision [] = {11, 8, 11};
  cpfloat_exponent_t emax [] = {15, 127, 127};
  size_t nformats = 3;
  printf("\n*** BINARY64 ***\n");
  const char outfile_conv_par[] = "./timing-clang-conv-par.dat";
  FILE *fid_conv_par = fopen(outfile_conv_par, "w");
  const char outfile_conv_seq[] = "./timing-clang-conv-seq.dat";
  FILE *fid_conv_seq = fopen(outfile_conv_seq, "w");
  const char outfile_conv_noalloc_seq[] = "./timing-clang-conv-noalloc-seq.dat";
  FILE *fid_conv_noalloc_seq = fopen(outfile_conv_noalloc_seq, "w");
  const char outfile_op_seq[] = "./timing-clang-op-seq.dat";
  FILE *fid_op_seq = fopen(outfile_op_seq, "w");
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    double *X = malloc(n * sizeof(*X));
    double *Y = malloc(n * sizeof(*Y));
    double *Xd, *Yd;
    for (j = 0; j < n; j++) { // generate normal numbers
      X[j] = fmin + rand() / (double)RAND_MAX;
      Y[j] = fmin + rand() / (double)RAND_MAX;
    }
    double *timing = malloc(NTESTS * sizeof(double));
    fprintf(fid_conv_par, "%6lu", sizes[i]);
    fprintf(fid_conv_seq, "%6lu", sizes[i]);
    fprintf(fid_conv_noalloc_seq, "%6lu", sizes[i]);
    fprintf(fid_op_seq, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (l = 0; l < nformats; l++) {
      fpopts->precision = precision[l];
      fpopts->emax = emax[l];

      /* Test parallel conversion with allocation. */
      for (k = 0; k < NTESTS; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        Xd = malloc(n * sizeof(*Xd));
        cpfloat(Xd, X, n, fpopts);
        clock_gettime(CLOCK_MONOTONIC, end);
        free(Xd);
        timing[k] = timedifference(start, end);
      }
      qsort(timing, NTESTS, sizeof(*timing), cmpfun);
      medtiming = timing[NTESTS / 2];
      fprintf(fid_conv_par, " %10.5e", medtiming);
      fprintf(stdout, " [P] %10.5e", medtiming);

      /* Test sequential conversion with allocation. */
      for (k = 0; k < NTESTS; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        Xd = malloc(n * sizeof(*Xd));
        cpfloat_sequential(Xd, Y, n, fpopts);
        clock_gettime(CLOCK_MONOTONIC, end);
        free(Xd);
        timing[k] = timedifference(start, end);
      }
      qsort(timing, NTESTS, sizeof(*timing), cmpfun);
      medtiming = timing[NTESTS/2];
      fprintf(fid_conv_seq, " %10.5e", medtiming);
      fprintf(stdout, " [C] %10.5e", medtiming);

      /* Test sequential conversion without allocation. */
      Xd = malloc(n * sizeof(*Xd));
      for (k = 0; k < NTESTS; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        cpfloat_sequential(Xd, Y, n, fpopts);
        clock_gettime(CLOCK_MONOTONIC, end);
        timing[k] = timedifference(start, end);
      }
      qsort(timing, NTESTS, sizeof(*timing), cmpfun);
      medtiming = timing[NTESTS/2];
      fprintf(fid_conv_noalloc_seq, " %10.5e", medtiming);
      fprintf(stdout, " [N] %10.5e", medtiming);

      /* Test sequential arithmetic operation. */
      Yd = malloc(n * sizeof(*Yd));
      double *Zd = malloc(n * sizeof(*Zd));
      for (k = 0; k < NTESTS; k++) {
        clock_gettime(CLOCK_MONOTONIC, start);
        cpf_add_sequential(Zd, Xd, Yd, n, fpopts);
        clock_gettime(CLOCK_MONOTONIC, end);
        timing[k] = timedifference(start, end);
      }
      free(Xd);
      free(Yd);
      free(Zd);
      qsort(timing, NTESTS, sizeof(*timing), cmpfun);
      medtiming = timing[NTESTS/2];
      fprintf(fid_op_seq, " %10.5e", medtiming);
      fprintf(stdout, " [A] %10.5e", medtiming);

      fprintf(stdout, " |");
    }

    free(X);
    free(Y);
    fprintf(fid_conv_par, "\n");
    fprintf(fid_conv_seq, "\n");
    fprintf(fid_conv_noalloc_seq, "\n");
    fprintf(fid_op_seq, "\n");
    fprintf(stdout, "\n");
  }

  fclose(fid_conv_par);
  fclose(fid_conv_seq);
  fclose(fid_conv_noalloc_seq);
  fclose(fid_op_seq);

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
