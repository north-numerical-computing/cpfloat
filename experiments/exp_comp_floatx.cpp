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

int main() {

  /* Input parameters. */
  size_t i, j, k, l, n;
  size_t ntests = 10;

  size_t nsizes = 28;
  size_t sizes [nsizes];
  size_t mult = 1;
  for (i = 1; i<=3; i++) {
    mult *= 10;
    for (j = 1; j<10; j++)
      sizes[9*(i-1)+(j-1)] = mult*j;
  }
  sizes[9*(i-1)] = mult*10;

  // size_t sizes [] = {10, 20, 50,
  //                    100, 200, 500,
  //                    1000, 2000, 5000,
  //                    10000};
  // size_t nsizes = 10;
  float fmin = ldexpf(1.,-14);

  struct timespec start, end;

  double * medtimings;
  medtimings = new double [nsizes];
  double * timing;
  timing = new double[ntests];

  /* FloatX */
  printf("\n*** FloatX compile-time ***\n");
  const char outfilefloatx [] = "./timing-floatx-seq.dat";
  FILE *fidx = fopen(outfilefloatx, "w");
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    double * Xd;
    Xd = new double [n];
    for (j=0; j<n; j++) // generate normal numbers
      Xd[j] = fmin + rand() / (double)RAND_MAX;
    fprintf(fidx, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);

    for (k=0; k<ntests; k++) {
      clock_gettime(CLOCK_MONOTONIC, &start);
      flx::floatx<5, 11, double> * Yd;
      Yd = new flx::floatx<5, 11, double> [n];
      for (j=0; j<n; j++)
        Yd[j] = Xd[j];
      clock_gettime(CLOCK_MONOTONIC, &end);
      delete[] (Yd);
      timing[k] = timedifference(&start, &end);
    }
    qsort(timing, ntests, sizeof(*timing), cmpfun);
    medtimings[i] = timing[ntests/2];
    fprintf(fidx, " %10.5e", medtimings[i]);
    fprintf(stdout, " %10.5e", medtimings[i]);

    for (k=0; k<ntests; k++) {
      clock_gettime(CLOCK_MONOTONIC, &start);
      flx::floatx<8, 8, double> * Yd;
      Yd = new flx::floatx<8, 8, double> [n];
      for (j=0; j<n; j++)
        Yd[j]= Xd[j];
      clock_gettime(CLOCK_MONOTONIC, &end);
      delete[] (Yd);
      timing[k] = timedifference(&start, &end);
    }
    qsort(timing, ntests, sizeof(*timing), cmpfun);
    medtimings[i] = timing[ntests/2];
    fprintf(fidx, " %10.5e", medtimings[i]);
    fprintf(stdout, " %10.5e", medtimings[i]);

    for (k=0; k<ntests; k++) {
      clock_gettime(CLOCK_MONOTONIC, &start);
      flx::floatx<8, 11,double> * Yd;
      Yd = new flx::floatx<8, 11, double> [n];
      for (j=0; j<n; j++)
        Yd[j]= Xd[j];
      clock_gettime(CLOCK_MONOTONIC, &end);
      delete[] (Yd);
      timing[k] = timedifference(&start, &end);
    }
    qsort(timing, ntests, sizeof(*timing), cmpfun);
    medtimings[i] = timing[ntests/2];
    fprintf(fidx, " %10.5e", medtimings[i]);
    fprintf(stdout, " %10.5e", medtimings[i]);

    fprintf(fidx, "\n");
    fprintf(stdout, "\n");
    delete[] (Xd);
  }
  fclose(fidx);

  /* floatxr */
  printf("\n*** FloatX runtime ***\n");
  size_t precision [] = {11, 8, 11};
  size_t exponent [] = {5, 8, 8};
  size_t nformats = 3;
  const char outfilefloatxr[] = "./timing-floatxr-seq.dat";
  FILE *fidxr = fopen(outfilefloatxr, "w");
  for (i = 0; i < nsizes; i++) {
    n = sizes[i] * sizes[i];
    double * Xd;
    Xd = new double [n];
    for (j=0; j<n; j++) // generate normal numbers
      Xd[j] = fmin + rand() / (double)RAND_MAX;
    fprintf(fidxr, "%6lu", sizes[i]);
    fprintf(stdout, "%6lu", sizes[i]);
    for (l=0; l<nformats; l++) {
      for (k=0; k<ntests; k++) {
        clock_gettime(CLOCK_MONOTONIC, &start);
        flx::floatxr<> ** Yd;
        Yd = new flx::floatxr<>* [n];
        for (j=0; j<n; j++)
          Yd[j]= new flx::floatxr<>(exponent[l], precision[l], Xd[j]);
        clock_gettime(CLOCK_MONOTONIC, &end);
        for (j=0; j<n; j++)
          delete (Yd[j]);
        delete[] (Yd);
        timing[k] = timedifference(&start, &end);
      }
      qsort(timing, ntests, sizeof(*timing), cmpfun);
      medtimings[i] = timing[ntests/2];
      fprintf(fidxr, " %10.5e", medtimings[i]);
      fprintf(stdout, " %10.5e", medtimings[i]);
    }
    fprintf(fidxr, "\n");
    fprintf(stdout, "\n");
    delete[] (Xd);
  }
         fclose(fidxr);

  delete[] (medtimings);
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
