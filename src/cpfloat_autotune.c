/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * This file is part of CPFloat.
 *
 * Running this program will update the threshold values in
 * cpfloat_threshold_binary32.h and cpfloat_threshold_binary64.h.
 */
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

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


bool parfaster_double(size_t n, optstruct *fpopts,
                      struct timespec *start, struct timespec *end,
                      double *seqtimings, double *partimings, size_t ntests) {
  size_t i;
  double *Xd = malloc(n * sizeof(*Xd));
  double *Yd = malloc(n * sizeof(*Yd));
  for (i = 0; i < n; i++)
    Xd[i] = rand() / (double)RAND_MAX;

  for (i = 0; i < ntests; i++) {
    clock_gettime(CLOCK_MONOTONIC, start);
    cpfloat_sequential(Yd, Xd, n, fpopts);
    clock_gettime(CLOCK_MONOTONIC, end);
    seqtimings[i] = timedifference(start, end);
    clock_gettime(CLOCK_MONOTONIC, start);
    cpfloat_parallel(Yd, Xd, n, fpopts);
    clock_gettime(CLOCK_MONOTONIC, end);
    partimings[i] = timedifference(start, end);
  }
  free(Xd);
  free(Yd);
  qsort(seqtimings, ntests, sizeof(*seqtimings), cmpfun);
  double seqtime = seqtimings[ntests/2];
  qsort(partimings, ntests, sizeof(*partimings), cmpfun);
  double partime = partimings[ntests/2];
  printf("[%7lu]   [%.2e, %.2e, %.2e]   [%.2e, %.2e, %.2e]\n",
         n, seqtimings[0], seqtime, seqtimings[ntests-1],
         partimings[0], partime, partimings[ntests-1]);
  return partime < seqtime ? true : false;
}

bool parfaster_float(size_t n, optstruct *fpopts,
                     struct timespec *start, struct timespec *end,
                     double *seqtimings, double *partimings, size_t ntests) {
  size_t i;
  float *Xd = malloc(n * sizeof(*Xd));
  float *Yd = malloc(n * sizeof(*Yd));
  for (i = 0; i < n; i++)
    Xd[i] = rand() / (float)RAND_MAX;
  for (i = 0; i < ntests; i++) {
    clock_gettime(CLOCK_MONOTONIC, start);
    cpfloatf_sequential(Yd, Xd, n, fpopts);
    clock_gettime(CLOCK_MONOTONIC, end);
    seqtimings[i] = timedifference(start, end);
    clock_gettime(CLOCK_MONOTONIC, start);
    cpfloatf_parallel(Yd, Xd, n, fpopts);
    clock_gettime(CLOCK_MONOTONIC, end);
    partimings[i] = timedifference(start, end);
  }
  free(Xd);
  free(Yd);
  qsort(seqtimings, ntests, sizeof(*seqtimings), cmpfun);
  float seqtime = seqtimings[ntests/2];
  qsort(partimings, ntests, sizeof(*partimings), cmpfun);
  float partime = partimings[ntests/2];
  printf("[%7lu]   [%.2e, %.2e, %.2e]   [%.2e, %.2e, %.2e]\n",
         n, seqtimings[0], seqtime, seqtimings[ntests-1],
         partimings[0], partime, partimings[ntests-1]);
  return partime < seqtime ? true : false;
}

int main() {

  /* Allocate fpopts struct and set fields to default. */
  static optstruct *fpopts;
  fpopts = malloc(sizeof(optstruct));
  strcpy(fpopts->format,"s");
  fpopts->precision = 24;
  fpopts->emax = 127;
  fpopts->subnormal = 0;
  fpopts->round = 1;
  fpopts->flip = 0;
  fpopts->p = 0.5;
  fpopts->explim = 1;

  size_t nmin, nmax, nmid;
  size_t ntests = 1000;
  struct timespec *start = malloc(sizeof(struct timespec));
  struct timespec *end = malloc(sizeof(struct timespec));

  double *seqtimings = malloc(ntests * sizeof(*seqtimings));
  double *partimings = malloc(ntests * sizeof(*partimings));

  int maxnumthreads = omp_get_max_threads();
  printf("Test using %d OpenMP threads.\n", maxnumthreads);

  char docstring [] =
    "/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */\n"
    "/* SPDX-License-Identifier: LGPL-2.1-or-later                         */\n"
    "\n/**\n"
    " * @file %s_threshold_%s.h\n"
    " * @brief Size of smallest `%s` array on which to use"
    " multiple OpenMP threads.\n"
    " */\n"
    "\n"
    "/**\n"
    " * @brief Size of smallest array on which %s() uses multiple threads.\n"
    " *\n"
    " * @details Threshold for switching between %s_sequential() and\n"
    " * %s_parallel() in %s(). The value of this constant is ignored\n"
    " * if the file that includes cpfloat_%s.h is compiled without OpenMP\n"
    " * support.\n"
    " */\n";

  /* Binary32 */
  nmin = 1;
  nmax = 1;
  while(!parfaster_float(nmax, fpopts, start, end,
                         seqtimings, partimings, ntests))
    nmax *= 2;
  nmid = (nmax + nmin) / 2;
  while(nmid != nmin && nmid != nmax) {
    if(parfaster_float(nmid, fpopts, start, end,
                       seqtimings, partimings, ntests))
      nmax = nmid;
    else
      nmin = nmid;
    nmid = (nmax + nmin) / 2;
  }
  printf("Optimal threshold for single is %zu.\n", nmax);
  const char filenamef [] = "./cpfloat_threshold_binary32.h";
  FILE *fidf = fopen(filenamef, "w");
  fprintf(fidf, docstring, "cpfloat", "binary32", "float",
          "cpfloatf", "cpfloatf", "cpfloatf", "cpfloatf", "binary32");
  fprintf(fidf, "#define OPENMP_THRESHOLD_float %zu", nmax);
  fclose(fidf);

  /* Binary64 */
  nmin = 1;
  nmax = 1;
  while(!parfaster_double(nmax, fpopts, start, end,
                          seqtimings, partimings, ntests))
    nmax *= 2;
  nmid = (nmax + nmin) / 2;
  while(nmid != nmin && nmid != nmax) {
    /* printf("[%5zu, %5zu, %5zu]\n", nmin, nmid, nmax); */
    if(parfaster_double(nmid, fpopts, start, end,
                        seqtimings, partimings, ntests))
      nmax = nmid;
    else
      nmin = nmid;
    nmid = (nmax + nmin) / 2;
  }
  printf("Optimal threshold for double is %zu.\n", nmax);
  const char filenamed [] = "./cpfloat_threshold_binary64.h";
  FILE *fidd = fopen(filenamed, "w");
  fprintf(fidd, docstring, "cpfloat", "binary64", "double",
          "cpfloat", "cpfloat", "cpfloat", "cpfloat", "binary64");
  fprintf(fidd, "#define OPENMP_THRESHOLD_double %zu", nmax);
  fclose(fidd);

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
