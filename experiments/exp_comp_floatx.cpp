/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * Compare new code with floatx.
 */

#include <iostream>

#include <cstdio>
#include <cstring>

#include "floatx.hpp"
#include "timing_fun.h"

template <short exp_bits, short sig_bits, typename backend_float>
void run_floatx_test(backend_float *Xd, size_t n, size_t ntests,
                     double *timing, FILE *fidx) {
  struct timespec start, end;
  for (size_t k = 0; k < ntests; k++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    flx::floatx<exp_bits, sig_bits, backend_float> *Yd;
    Yd = new flx::floatx<exp_bits, sig_bits, backend_float>[n];
    for (size_t j = 0; j < n; j++)
      Yd[j] = Xd[j];
    clock_gettime(CLOCK_MONOTONIC, &end);
    delete[](Yd);
    timing[k] = timedifference(&start, &end);
  }
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  double medtiming = timing[ntests / 2];
  fprintf(fidx, " %10.5e", medtiming);
  fprintf(stdout, " %10.5e", medtiming);
}

template <typename backend_float>
void run_floatxr_test(short exp_bits, short sig_bits,
                      backend_float *Xd, size_t n, size_t ntests,
                      double *timing, FILE *fidx) {
  struct timespec start, end;
  for (size_t k = 0; k < ntests; k++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    flx::floatxr<double, short> **Yd;
    Yd = new flx::floatxr<double, short> *[n];
    for (size_t j = 0; j < n; j++)
      Yd[j] = new flx::floatxr<double, short>(exp_bits, sig_bits, Xd[j]);
    clock_gettime(CLOCK_MONOTONIC, &end);
    for (size_t j = 0; j < n; j++)
      delete (Yd[j]);
    delete[](Yd);
    timing[k] = timedifference(&start, &end);
  }
  qsort(timing, ntests, sizeof(*timing), cmpfun);
  double medtiming = timing[ntests / 2];
  fprintf(fidx, " %10.5e", medtiming);
  fprintf(stdout, " %10.5e", medtiming);
}

int main() {

  /* Input parameters. */
  const size_t ntests = 10;
  const size_t nsizes = 28;
  size_t sizes [nsizes];
  size_t mult = 1;
  size_t i, j;
  for (i = 1; i<=3; i++) {
    mult *= 10;
    for (j = 1; j<10; j++)
      sizes[9*(i-1)+(j-1)] = mult*j;
  }
  sizes[9*(i-1)] = mult*10;

  constexpr size_t precision[3] = {11, 8, 11};
  constexpr size_t exponent[3] = {5, 8, 8};
  constexpr size_t nformats = 3;

  // Smallest normal in binary16. Adding this to numbers generated in [0,1]
  // ensures that values are normal in binary16, bfloat16, and TensorFloat-32.
  double fmin = ldexpf(1.,-14);

  double *timing;
  timing = new double[ntests];

  /* Compile-time FloatX values (floatx class). */
  printf("\n*** FloatX compile-time ***\n");
  const char outfilefloatx [] = "./timing-floatx-seq.dat";
  FILE *fidx = fopen(outfilefloatx, "w");
  for (i = 0; i < nsizes; i++) {
    size_t n = sizes[i] * sizes[i];
    double * Xd;
    Xd = new double [n];
    for (j=0; j<n; j++) // generate normal numbers
      Xd[j] = fmin + rand() / (double)RAND_MAX;
    fprintf(fidx, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);

    // binary16
    run_floatx_test<5, 11, double>(Xd, n, ntests, timing, fidx);
    // bfloat16
    run_floatx_test<8, 8, double>(Xd, n, ntests, timing, fidx);
    // TensorFloat-32
    run_floatx_test<8, 11, double>(Xd, n, ntests, timing, fidx);

    fprintf(fidx, "\n");
    fprintf(stdout, "\n");
    delete[](Xd);
  }
  fclose(fidx);

  /* Run-time FloatX values (floatxr class). */
  printf("\n*** FloatX runtime ***\n");
  const char outfilefloatxr[] = "./timing-floatxr-seq.dat";
  FILE *fidxr = fopen(outfilefloatxr, "w");
  for (i = 0; i < nsizes; i++) {
    size_t n = sizes[i] * sizes[i];
    double * Xd;
    Xd = new double [n];
    for (j=0; j<n; j++) // generate normal numbers
      Xd[j] = fmin + rand() / (double)RAND_MAX;
    fprintf(fidxr, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (size_t l=0; l<nformats; l++) {
      run_floatxr_test<double>(exponent[l], precision[l],
                               Xd, n, ntests, timing, fidx);
    }
    fprintf(fidxr, "\n");
    fprintf(stdout, "\n");
    delete[] (Xd);
  }
         fclose(fidxr);

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
