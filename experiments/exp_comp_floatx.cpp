/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * Compare CPFloat with floatx.
 */

#include <iostream>

#include <cstdio>
#include <cstring>

#include "floatx.hpp"
#include "timing_fun.h"

template <short exp_bits, short sig_bits, typename backend_float>
void run_floatx_test(backend_float *X, backend_float *Y, size_t n,
                     size_t ntests, double *timing,
                     FILE *fidx_conv, FILE *fidx_conv_noalloc, FILE *fidx_op) {
  struct timespec start, end;
  double medtiming;

  // Test conversion with allocation.
  for (size_t k = 0; k < ntests; k++) {
    // flx::floatx<exp_bits, sig_bits, backend_float> *Yd;
    clock_gettime(CLOCK_MONOTONIC, &start);
    auto Xd = new flx::floatx<exp_bits, sig_bits, backend_float>[n];
    for (size_t j = 0; j < n; j++)
      Xd[j] = X[j];
    clock_gettime(CLOCK_MONOTONIC, &end);
    delete[](Xd);
    timing[k] = timedifference(&start, &end);
  }
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  medtiming = timing[ntests / 2];
  fprintf(fidx_conv, " %10.5e", medtiming);
  fprintf(stdout, " [C] %10.5e", medtiming);

  // Test conversion without allocation.
  auto Xd = new flx::floatx<exp_bits, sig_bits, backend_float>[n];
  for (size_t k = 0; k < ntests; k++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (size_t j = 0; j < n; j++)
      Xd[j] = X[j];
    clock_gettime(CLOCK_MONOTONIC, &end);
    timing[k] = timedifference(&start, &end);
  }
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  medtiming = timing[ntests / 2];
  fprintf(fidx_conv_noalloc, " %10.5e", medtiming);
  fprintf(stdout, " [N] %10.5e", medtiming);

  // Test arithmetic operation.
  auto Yd = new flx::floatx<exp_bits, sig_bits, backend_float>[n];
  auto Zd = new flx::floatx<exp_bits, sig_bits, backend_float>[n];
  for (size_t j = 0; j < n; j++)
    Yd[j] = Y[j];
  for (size_t k = 0; k < ntests; k++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (size_t j = 0; j < n; j++)
      Zd[j] = Yd[j] + Zd[j];
    clock_gettime(CLOCK_MONOTONIC, &end);
    timing[k] = timedifference(&start, &end);
  }
  delete[](Xd);
  delete[](Yd);
  delete[](Zd);
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  medtiming = timing[ntests / 2];
  fprintf(fidx_op, " %10.5e", medtiming);
  fprintf(stdout, " [A] %10.5e", medtiming);
}

template <typename backend_float>
void run_floatxr_test(short exp_bits, short sig_bits,
                      backend_float *X, backend_float *Y, size_t n,
                      size_t ntests, double *timing,
                      FILE *fidx_conv, FILE *fidx_conv_noalloc, FILE *fidx_op) {
  struct timespec start, end;
  double medtiming;

  // Test conversion with allocation.
  for (size_t k = 0; k < ntests; k++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    auto Xd = new flx::floatxr<double, short> *[n];
    for (size_t j = 0; j < n; j++)
      Xd[j] = new flx::floatxr<double, short>(exp_bits, sig_bits);
    for (size_t j = 0; j < n; j++)
      *Xd[j] = X[j];
    clock_gettime(CLOCK_MONOTONIC, &end);
    for (size_t j = 0; j < n; j++)
      delete (Xd[j]);
    delete[](Xd);
    timing[k] = timedifference(&start, &end);
  }
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  medtiming = timing[ntests / 2];
  fprintf(fidx_conv, " %10.5e", medtiming);
  fprintf(stdout, " [C] %10.5e", medtiming);

  // Test conversion without allocation.
  auto Xd = new flx::floatxr<double, short> *[n];
  for (size_t j = 0; j < n; j++)
    Xd[j] = new flx::floatxr<double, short>(exp_bits, sig_bits);
  for (size_t k = 0; k < ntests; k++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (size_t j = 0; j < n; j++)
      *Xd[j] = X[j];
    clock_gettime(CLOCK_MONOTONIC, &end);
    timing[k] = timedifference(&start, &end);
  }
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  medtiming = timing[ntests / 2];
  fprintf(fidx_conv_noalloc, " %10.5e", medtiming);
  fprintf(stdout, " [C] %10.5e", medtiming);

  // Test arithmetic operation.
  auto Yd = new flx::floatxr<double, short> *[n];
  auto Zd = new flx::floatxr<double, short> *[n];
  for (size_t j = 0; j < n; j++) {
    Yd[j] = new flx::floatxr<double, short>(exp_bits, sig_bits, Y[j]);
    Zd[j] = new flx::floatxr<double, short>(exp_bits, sig_bits);
  }
  for (size_t k = 0; k < ntests; k++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (size_t j = 0; j < n; j++)
      *Zd[j] = *Xd[j] + *Yd[j];
    clock_gettime(CLOCK_MONOTONIC, &end);
    timing[k] = timedifference(&start, &end);
  }
    for (size_t j = 0; j < n; j++) {
    delete (Xd[j]);
    delete (Yd[j]);
    delete (Zd[j]);
  }
  delete[](Xd);
  delete[](Yd);
  delete[](Zd);
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  medtiming = timing[ntests / 2];
  fprintf(fidx_op, " %10.5e", medtiming);
  fprintf(stdout, " [A] %10.5e", medtiming);
}





int main() {

  /* Input parameters. */
  const size_t nsizes = 28;
  size_t sizes[nsizes];
  size_t mult = 1;
  size_t i, j;
  for (i = 1; i <= 3; i++) {
    mult *= 10;
    for (j = 1; j < 10; j++)
      sizes[9 * (i - 1) + (j - 1)] = mult * j;
  }
  sizes[9 * (i - 1)] = mult * 10;

  constexpr size_t precision[3] = {11, 8, 11};
  constexpr size_t exponent[3] = {5, 8, 8};
  constexpr size_t nformats = 3;

  // Smallest normal in binary16. Adding this to numbers generated in [0,1]
  // ensures that values are normal in binary16, bfloat16, and TensorFloat-32.
  double fmin = ldexpf(1., -14);

  double *timing;
  timing = new double[NTESTS];

  /* Compile-time FloatX values (floatx class). */
  printf("\n*** FloatX compile-time ***\n");
  const char outfile_fx_conv[] = "./timing-floatx-conv-seq.dat";
  FILE *fidx_conv = fopen(outfile_fx_conv, "w");
  const char outfile_fx_conv_noalloc[] = "./timing-floatx-conv-noalloc-seq.dat";
  FILE *fidx_conv_noalloc = fopen(outfile_fx_conv_noalloc, "w");
  const char outfile_fx_op[] = "./timing-floatx-op-seq.dat";
  FILE *fidx_op = fopen(outfile_fx_op, "w");
  for (i = 0; i < nsizes; i++) {
    size_t n = sizes[i] * sizes[i];
    double *X = new double[n];
    double *Y = new double[n];
    for (j = 0; j < n; j++) { // generate normal numbers
      X[j] = fmin + rand() / (double)RAND_MAX;
      Y[j] = fmin + rand() / (double)RAND_MAX;
    }
    fprintf(fidx_conv, "%6lu", sizes[i]);
    fprintf(fidx_conv_noalloc, "%6lu", sizes[i]);
    fprintf(fidx_op, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);

    // binary16
    run_floatx_test<5, 11, double>(X, Y, n, NTESTS, timing,
                                   fidx_conv, fidx_conv_noalloc, fidx_op);
    fprintf(stdout, " |");
    // bfloat16
    run_floatx_test<8, 8, double>(X, Y, n, NTESTS, timing,
                                  fidx_conv, fidx_conv_noalloc, fidx_op);
    fprintf(stdout, " |");
    // TensorFloat-32
    run_floatx_test<8, 11, double>(X, Y, n, NTESTS, timing,
                                   fidx_conv, fidx_conv_noalloc, fidx_op);
    fprintf(stdout, " |");

    fprintf(fidx_conv, "\n");
    fprintf(fidx_conv_noalloc, "\n");
    fprintf(fidx_op, "\n");
    fprintf(stdout, "\n");
    delete[](X);
    delete[](Y);
  }
  fclose(fidx_conv);
  fclose(fidx_conv_noalloc);
  fclose(fidx_op);

  /* Run-time FloatX values (floatxr class). */
  printf("\n*** FloatX runtime ***\n");
  const char outfile_fxr_conv [] = "./timing-floatxr-conv-seq.dat";
  FILE *fidxr_conv = fopen(outfile_fxr_conv, "w");
  const char outfile_fxr_conv_noalloc[] = "./timing-floatxr-conv-noalloc-seq.dat";
  FILE *fidxr_conv_noalloc = fopen(outfile_fxr_conv_noalloc, "w");
  const char outfile_fxr_op[] = "./timing-floatxr-op-seq.dat";
  FILE *fidxr_op = fopen(outfile_fxr_op, "w");
  for (i = 0; i < nsizes; i++) {
    size_t n = sizes[i] * sizes[i];
    double *X = new double[n];
    double *Y = new double[n];
    for (j = 0; j < n; j++) { // generate normal numbers
      X[j] = fmin + rand() / (double)RAND_MAX;
      Y[j] = fmin + rand() / (double)RAND_MAX;
    }
    fprintf(fidxr_conv, "%6lu", sizes[i]);
    fprintf(fidxr_conv_noalloc, "%6lu", sizes[i]);
    fprintf(fidxr_op, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (size_t l=0; l<nformats; l++) {
      run_floatxr_test<double>(exponent[l], precision[l],
                               X, Y, n, NTESTS, timing,
                               fidxr_conv, fidxr_conv_noalloc, fidxr_op);
      fprintf(stdout, " |");
    }
    fprintf(fidxr_conv, "\n");
    fprintf(fidxr_conv_noalloc, "\n");
    fprintf(fidxr_op, "\n");
    fprintf(stdout, "\n");
    delete[](X);
    delete[](Y);
  }
  fclose(fidxr_conv);
  fclose(fidxr_conv_noalloc);
  fclose(fidxr_op);

  delete[] (timing);

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
