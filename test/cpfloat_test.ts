/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include "cpfloat_binary32.h"
#include "cpfloat_binary64.h"

#define MINFORMAT 0
#define MAXFORMAT 2
#define MINMODE  -1
#define MAXMODE   8
#define NREPS  1000

/* Define default target formats. */
// binary16, bfloat16, TensorFloat-32
static size_t precision [] = {11, 8, 11};
static size_t emax [] = {15, 127, 127};
static size_t nformats = 3;

/* Structure for options and fixtures. */
optstruct *fpopts;

void fpopts_setup(void) {
  fpopts = malloc(sizeof(optstruct));
  fpopts->flip = 0;
  fpopts->p = 0.5;
  fpopts->bitseed = NULL;
  fpopts->randseedf = NULL;
  fpopts->randseed = NULL;
}

void fpopts_teardown(void) {
  free(fpopts);
}

#define INTOFd(x)(((fpint){.fpval = (double)(x)}).intval)
#define FPOFd(x)(((fpint){.intval = (uint64_t)(x)}).fpval)

#define INTOFf(x)(((fpintf){.fpval = (float)(x)}).intval)
#define FPOFf(x)(((fpintf){.intval = (uint32_t)(x)}).fpval)

/***********************
 * Auxiliary functions *
 ***********************/

/* Return values of interest */
static inline
uint64_t intminnormal_double(optstruct *fpopts) {
  return INTOFd(ldexp(1., 1-fpopts->emax));
}

static inline
uint32_t intminnormal_float(optstruct *fpopts) {
  return INTOFf(ldexp(1., 1-fpopts->emax));
}

static inline
double minsubnormal(optstruct *fpopts) {
  return ldexp(1., 2 - fpopts->emax - fpopts->precision);
}

static inline
double maxsubnormal(optstruct *fpopts) {
  return ldexp(1., 1 - fpopts->emax) *
    (1 - ldexp(1., 1 - fpopts->precision));
}

static inline
double minnormal(optstruct *fpopts) {
  return ldexp(1., 1-fpopts->emax);
}

static inline
double maxnormal(optstruct *fpopts) {
  return ldexp(1., fpopts->emax) *
    (2 - ldexp(1., 1 - fpopts->precision));
}

static inline
double maxbound(optstruct *fpopts) {
  return ldexp(1., fpopts->emax) *
    (2 - ldexp(1., -fpopts->precision));
}


static inline
double inf_double() {
  return FPOFd(0x7ff0000000000000);
}

static inline
float inf_float() {
  return FPOFf(0x7f800000);
}

static inline
double nan_double() {
  return FPOFd(0x7ff8000000000000);
}

static inline
float nan_float() {
  return FPOFf(0x7fa00000);
}

static inline
void init_intarray_double(uint64_t *x, size_t n,
                          uint64_t first, uint64_t step) {
  x[0] = first;
  for (size_t i = 1; i < n; i++)
    x[i] = x[i-1] + step;
}

static inline
void init_intarray_float(uint32_t *x, size_t n,
                         uint32_t first, uint32_t step) {
  x[0] = first;
  for (size_t i = 1; i < n; i++)
    x[i] = x[i-1] + step;
}

static inline
void init_intarray_rounding_double(uint64_t *x, size_t n,
                                   uint64_t first, uint64_t step) {
  x[0] = first + 1ul;
  x[1] = first + (step >> 1);
  for (size_t i = 2; i < n; i+=3) {
    first += step;
    x[i] = first - 1ul;
    if (i+1 < n)
      x[i+1] = first + 1ul;
    if (i+2 < n)
      x[i+2] = first + (step >> 1);
  }
}

static inline
void init_intarray_rounding_float(uint32_t *x, size_t n,
                                  uint32_t first, uint32_t step) {
  x[0] = first + 1u;
  x[1] = first + (step >> 1);
  for (size_t i = 2; i < n; i+=3) {
    first += step;
    x[i] = first - 1u;
    if (i+1 < n)
      x[i+1] = first + 1ul;
    if (i+2 < n)
      x[i+2] = first + (step >> 1);
  }
}

static inline
void init_fparray_double(double *x, size_t n,
                         double first, double step) {
  x[0] = first;
  for (size_t i = 1; i < n; i++)
    x[i] = x[i-1] + step;
}

static inline
void init_fparray_float(float *x, size_t n,
                        float first, float step) {
  x[0] = first;
  for (size_t i = 1; i < n; i++)
    x[i] = x[i-1] + step;
}

static inline
void init_fparray_rounding_double(double *x, size_t n,
                                  double first, double step) {
  x[0] = FPOFd(INTOFd(first) + 1lu);
  x[1] = first + step/2;
  for (size_t i = 2; i < n; i+=3) {
    first += step;
    x[i] = FPOFd(INTOFd(first) - 1lu);
    if (i+1 < n)
        x[i+1] = FPOFd(INTOFd(first) + 1lu);
    if (i+2 < n)
        x[i+2] = first + step/2;
  }
}

static inline
void init_fparray_rounding_float(float *x, size_t n,
                                 float first, float step) {
  x[0] = FPOFf(INTOFf(first) + 1lu);
  x[1] = first + step/2;
  for (size_t i = 2; i < n; i+=3) {
    first += step;
    x[i] = FPOFf(INTOFf(first) - 1lu);
    if (i+1 < n)
        x[i+1] = FPOFf(INTOFf(first) + 1lu);
    if (i+2 < n)
        x[i+2] = first + step/2;
  }
}

static inline
void init_intarray_rounding_stoc_double(uint64_t *x, size_t n,
                                        uint64_t first, uint64_t step) {
  x[0] = first + (step >> 2);
  x[1] = x[0] + (step >> 2);
  x[2] = x[1] + (step >> 2);
  for (size_t i = 3; i < n; i++) {
    x[i] = x[i-3] + step;
  }
}

static inline
void init_intarray_rounding_stoc_float(uint32_t *x, size_t n,
                                       uint32_t first, uint32_t step) {
  x[0] = first + (step >> 2);
  x[1] = x[0] + (step >> 2);
  x[2] = x[1] + (step >> 2);
  for (size_t i = 3; i < n; i++) {
    x[i] = x[i-3] + step;
  }
}


static inline
void init_fparray_rounding_stoc_double(double *x, size_t n,
                                       double first, double step) {
  x[0] = first + step/4;
  x[1] = x[0] + step/4;
  x[2] = x[1] + step/4;
  for (size_t i = 3; i < n; i++) {
    x[i] = x[i-3] + step;
  }
}

static inline
void init_fparray_rounding_stoc_float(float *x, size_t n,
                                      float first, float step) {
  x[0] = first + step/4;
  x[1] = x[0] + step/4;
  x[2] = x[1] + step/4;
  for (size_t i = 3; i < n; i++) {
    x[i] = x[i-3] + step;
  }
}

static inline
void csign_intarray_double(uint64_t *x, size_t n) {
  for (size_t i = 0; i < n; i++)
    x[i] = x[i] ^ 0x8000000000000000;
}

static inline
void csign_intarray_float(uint32_t *x, size_t n) {
  for (size_t i = 0; i < n; i++)
    x[i] = x[i] ^ 0x80000000;
}

static inline
void invprob_array_double(double *x, size_t n) {
  for (size_t i = 0; i < n; i++)
    x[i] = 1 - x[i];
}

/* Comparison with DWIM for NaN values. */
static inline
int nan_safe_compare_double(double a, double b) {
  return a == b || (isnan(a) && isnan(b));
}

static inline
int nan_safe_compare_float(float a, float b) {
  return a == b || (isnan(a) && isnan(b));
}

/* Run on single array and check results. */
static inline
void check_array_double(double *y, double *x, double *ref,
                        size_t n, optstruct * fpopts) {
  cpfloat(y, x, n, fpopts);
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_double(y[j], ref[j])) {
      printf("DOUBLE\n");
      printf("***\nj = %ld\nin  = %23.15e [%lX]\nout = %23.15e [%lX]\nref = %23.15e [%lX]\nr = %d\npr = %d\ne = %d\ns = %d\n",
             j,
             x[j], *(uint64_t *)(x+j),
             y[j], *(uint64_t *)(y+j),
             ref[j], *(uint64_t *)(ref+j),
             fpopts->round, fpopts->precision, fpopts->emax,
             fpopts->subnormal);
    }
    ck_assert(nan_safe_compare_double(y[j], ref[j]));
  }
}

static inline
void check_array_float(float *y, float *x, float *ref,
                       size_t n, optstruct * fpopts) {
  cpfloatf(y, x, n, fpopts);
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_float(y[j], ref[j])) {
      printf("FLOAT\n");
      printf(
          "***\nj = %ld\nin  = %16.8e [%X]\nout = %16.8e [%X]\nref = %16.8e [%X]\nr = %d\npr = %d\ne = %d\ns = %d\n",
          j,
          x[j], * (uint32_t *)(x+j),
          y[j], * (uint32_t *)(y+j),
          ref[j], * (uint32_t *)(ref+j),
          fpopts->round, fpopts->precision,
          fpopts->emax, fpopts->subnormal);
    }
    ck_assert(nan_safe_compare_float(y[j], ref[j]));
  }
}

static inline
void check_array_stoc_double(double *tmpin, double *tmpout,
                             double *x, double *prounddown,
                             size_t n, optstruct *fpopts) {
  char mode = fpopts->round;
  for (size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < NREPS; j++)
      tmpin[j] = x[i];
    cpfloat(tmpout, tmpin, NREPS, fpopts);
    double counter [] = {0, 0};
    double xup = 0;
    double xdown = 0;
    fpopts->round = 2;
    cpfloat(&xup, x+i, 1, fpopts);
    fpopts->round = 3;
    cpfloat(&xdown, x+i, 1, fpopts);
    for (size_t j = 0; j < NREPS; j++)
      if (tmpout[j] == xup)
        counter[1]++;
      else if (tmpout[j] == xdown)
        counter[0]++;
      else
        ck_abort_msg("Not rounding to either closest number.");
    fpopts->round = mode;
    if (fabs(counter[0]/NREPS - prounddown[i % 3]) > 0.1) {
      ck_abort_msg("Error in stochasting rounding.");
    }
  }
}

static inline
void check_array_stoc_float(float *tmpin, float *tmpout,
                            float *x, double *prounddown,
                            size_t n, optstruct *fpopts) {
  char mode = fpopts->round;
  for (size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < NREPS; j++)
      tmpin[j] = x[i];
    cpfloatf(tmpout, tmpin, NREPS, fpopts);
    double counter [] = {0, 0};
    float xup = 0;
    float xdown = 0;
    fpopts->round = 2;
    cpfloatf(&xup, x+i, 1, fpopts);
    fpopts->round = 3;
    cpfloatf(&xdown, x+i, 1, fpopts);
    for (size_t j = 0; j < NREPS; j++)
      if (tmpout[j] == xup)
        counter[1]++;
      else if (tmpout[j] == xdown)
        counter[0]++;
      else
        ck_abort_msg("Not rounding to either closest number.");
    fpopts->round = mode;
    if (fabs(counter[0]/(double)NREPS - prounddown[i % 3]) > 0.1) {
      ck_abort_msg("Error in stochasting rounding.");
    }
  }
}

static inline
void check_array_equi_double(double *tmpin, double *tmpout,
                             double *x, double *prounddown,
                             size_t n, optstruct *fpopts) {
  for (size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < NREPS; j++)
      tmpin[j] = x[i];
    cpfloat(tmpout, tmpin, NREPS, fpopts);
    double counter [] = {0, 0};
    for (size_t j = 0; j < NREPS; j++)
      if (tmpout[j] > x[i])
        counter[1]++;
      else
        counter[0]++;
    if (fabs(counter[0]/NREPS - *prounddown) > 0.1) {
      ck_abort_msg("Error in stochasting rounding.");
    }
  }
}

static inline
void check_array_equi_float(float *tmpin, float *tmpout,
                            float *x, double *prounddown,
                            size_t n, optstruct *fpopts) {
  for (size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < NREPS; j++)
      tmpin[j] = x[i];
    cpfloatf(tmpout, tmpin, NREPS, fpopts);
    double counter [] = {0, 0};
    for (size_t j = 0; j < NREPS; j++)
      if (tmpout[j] > x[i])
        counter[1]++;
      else
        counter[0]++;
    if (fabs(counter[0]/(double)NREPS - *prounddown) > 0.1) {
      ck_abort_msg("Error in stochasting rounding.");
    }
  }
}

void select_tests_det_double(double *y, double *x, double *ref,
                             size_t n, optstruct * fpopts,
                             int minmode, int maxmode,
                             size_t minformat, size_t maxformat,
                             int subnormalin, int explimin) {

  bool allsubnormal = (subnormalin != 0 && subnormalin != 1);
  size_t minsubnormal = allsubnormal ? 0 : subnormalin;
  size_t maxsubnormal = allsubnormal ? 1 : subnormalin;

  bool allexplim = (explimin != 0 && explimin != 1);
  size_t minexplim = allexplim ? 0 : explimin;
  size_t maxexplim = allexplim ? 1 : explimin;

  int round;
  size_t format, subnormal, explim;
  for (round = minmode; round <= maxmode; round++) {
    fpopts->round = round;
    for (format = minformat; format <= maxformat; format++) {
      fpopts->precision = precision[format];
      fpopts->emax = emax[format];
      for (subnormal = minsubnormal; subnormal <= maxsubnormal; subnormal++) {
        fpopts->subnormal = subnormal;
        for (explim = minexplim; explim <= maxexplim; explim++) {
          fpopts->explim = explim;
          check_array_double(y, x, ref, n, fpopts);
        }
      }
    }
  }
}

void select_tests_det_float(float *y, float *x, float *ref,
                            size_t n, optstruct * fpopts,
                            int minmode, int maxmode,
                            size_t minformat, size_t maxformat,
                            int subnormalin, int explimin) {

  bool allsubnormal = (subnormalin != 0 && subnormalin != 1);
  size_t minsubnormal = allsubnormal ? 0 : subnormalin;
  size_t maxsubnormal = allsubnormal ? 1 : subnormalin;

  bool allexplim = (explimin != 0 && explimin != 1);
  size_t minexplim = allexplim ? 0 : explimin;
  size_t maxexplim = allexplim ? 1 : explimin;

  int round;
  size_t format, subnormal, explim;
  for (round = minmode; round <= maxmode; round++) {
    fpopts->round = round;
    for (format = minformat; format <= maxformat; format++) {
      fpopts->precision = precision[format];
      fpopts->emax = emax[format];
      for (subnormal = minsubnormal; subnormal <= maxsubnormal; subnormal++) {
        fpopts->subnormal = subnormal;
        for (explim = minexplim; explim <= maxexplim; explim++) {
          fpopts->explim = explim;
          check_array_float(y, x, ref, n, fpopts);
        }
      }
    }
  }
}

void select_tests_stoc_double(double *tmpin, double *tmpout,
                              double *x, double *prounddown,
                              size_t n, optstruct * fpopts,
                              int mode,
                              size_t minformat, size_t maxformat,
                              int subnormalin, int explimin,
                              size_t testmode) {

  bool allsubnormal = (subnormalin != 0 && subnormalin != 1);
  size_t minsubnormal = allsubnormal ? 0 : subnormalin;
  size_t maxsubnormal = allsubnormal ? 1 : subnormalin;

  bool allexplim = (explimin != 0 && explimin != 1);
  size_t minexplim = allexplim ? 0 : explimin;
  size_t maxexplim = allexplim ? 1 : explimin;

  size_t format, subnormal, explim;
  fpopts->round = mode;
  for (format = minformat; format <= maxformat; format++) {
    fpopts->precision = precision[format];
    fpopts->emax = emax[format];
    for (subnormal = minsubnormal; subnormal <= maxsubnormal; subnormal++) {
      fpopts->subnormal = subnormal;
      for (explim = minexplim; explim <= maxexplim; explim++) {
        fpopts->explim = explim;
        double *outarray = testmode == 1? tmpout : tmpin;
        if (mode == 5)
          check_array_stoc_double(tmpin, outarray, x, prounddown, n, fpopts);
        else if (mode == 6)
          check_array_equi_double(tmpin, outarray, x, prounddown, n, fpopts);
      }
    }
  }
}

void select_tests_stoc_float(float *tmpin, float *tmpout,
                             float *x, double *prounddown,
                             size_t n, optstruct * fpopts,
                             int mode,
                             size_t minformat, size_t maxformat,
                             int subnormalin, int explimin,
                             size_t testmode) {

  bool allsubnormal = (subnormalin != 0 && subnormalin != 1);
  size_t minsubnormal = allsubnormal ? 0 : subnormalin;
  size_t maxsubnormal = allsubnormal ? 1 : subnormalin;

  bool allexplim = (explimin != 0 && explimin != 1);
  size_t minexplim = allexplim ? 0 : explimin;
  size_t maxexplim = allexplim ? 1 : explimin;

  size_t format, subnormal, explim;
  fpopts->round = mode;
  for (format = minformat; format <= maxformat; format++) {
    fpopts->precision = precision[format];
    fpopts->emax = emax[format];
    for (subnormal = minsubnormal; subnormal <= maxsubnormal; subnormal++) {
      fpopts->subnormal = subnormal;
      for (explim = minexplim; explim <= maxexplim; explim++) {
        fpopts->explim = explim;
        float *outarray = testmode == 1? tmpout : tmpin;
        if (mode == 5)
          check_array_stoc_float(tmpin, outarray, x, prounddown, n, fpopts);
        else if (mode == 6)
          check_array_equi_float(tmpin, outarray, x, prounddown, n, fpopts);
      }
    }
  }
}

double *alloc_init_array_double(double *x, size_t n) {
  double *y = malloc(n * sizeof(*y));
  for (size_t i=0; i<n; i++)
    y[i] = x[i];
  return y;
}

float *alloc_init_array_float(float *x, size_t n) {
  float *y = malloc(n * sizeof(*y));
  for (size_t i=0; i<n; i++)
    y[i] = x[i];
  return y;
}

void copy_array_double(double *y, double *x, size_t n) {
  for (size_t i=0; i<n; i++)
    y[i] = x[i];
}

void copy_array_float(float *y, float *x, size_t n) {
  for (size_t i=0; i<n; i++)
    y[i] = x[i];
}

double *allocate_array_double(double *x, size_t n, size_t mode) {
  double *y;
  if (mode == 1)
    y = malloc(n * sizeof(*y));
  else
    y = x;
  return y;
}

float *allocate_array_float(float *x, size_t n, size_t mode) {
  float *y;
  if (mode == 1)
    y = malloc(n * sizeof(*y));
  else
    y = x;
  return y;
}

void free_array_double(double *y, size_t mode) {
  if (mode == 1)
    free(y);
}

void free_array_float(float *y, size_t mode) {
  if (mode == 1)
    free(y);
}

void reset_array_double(double *y, size_t n, size_t mode) {
  if (mode == 1)
    for (size_t i=0; i<n; i++)
      y[i] = NAN;
}

void reset_array_float(float *y, size_t n, size_t mode) {
  if (mode == 1)
    for (size_t i=0; i<n; i++)
      y[i] = NAN;
}



/*********
 * TESTS *
 *********/
#suite cpfloat

/* 1. Values that do not need rounding. */
#tcase no_rounding

#test no_rounding_special_values
printf("1a. No rounding: special values\n");
for (size_t mode=1; mode<3; mode++) {
  /* Infinity and non-a-number. */
  size_t n = 4;
  double xd [] = {0, -0, inf_double(), -inf_double(), nan_double()};
  double *zd = alloc_init_array_double(xd, n);
  double *yd = allocate_array_double(xd, n, mode);
  select_tests_det_double(yd, xd, zd, n, fpopts,
                          MINMODE, MAXMODE, 0, 2, -1, -1);
  free(zd);
  free_array_double(yd, mode);

  float xf [] = {0, -0, inf_float(), -inf_float(), nan_float()};
  float *zf = alloc_init_array_float(xf, n);
  float *yf = allocate_array_float(xf, n, mode);
  select_tests_det_float(yf, xf, zf, n, fpopts,
                         MINMODE, MAXMODE, 0, 2, -1, -1);
  free(zf);
  free_array_float(yf, mode);
}

/* Exactly representable subnormals. */
#test no_rounding_subnormal_numbers
printf("1b. No rounding: subnormal numbers\n");
for (size_t mode=1; mode<3; mode++) {
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    size_t n = ldexp(1.,fpopts->precision-1) - 1; // number of subnormals
    double *xd = malloc(n * sizeof(*xd));
    init_fparray_double(xd, n, minsubnormal(fpopts), minsubnormal(fpopts));
    double *zd = alloc_init_array_double(xd, n);
    double *yd = allocate_array_double(xd, n, mode);
    select_tests_det_double(yd, xd, zd, n, fpopts, MINMODE, MAXMODE,
                            i, i, 1, 0);
    csign_intarray_double((uint64_t *)xd, n);
    csign_intarray_double((uint64_t *)zd, n);
    select_tests_det_double(yd, xd, zd, n, fpopts, MINMODE, MAXMODE,
                            i, i, 1, 0);
    free(xd);
    free(zd);
    free_array_double(yd, mode);

    float *xf = malloc(n * sizeof(*xf));
    init_fparray_float(xf, n, minsubnormal(fpopts), minsubnormal(fpopts));
    float *zf = alloc_init_array_float(xf, n);
    float *yf = allocate_array_float(xf, n, mode);
    select_tests_det_float(yf, xf, zf, n, fpopts, MINMODE, MAXMODE,
                           i, i, 1, 0);
    csign_intarray_float((uint32_t *)xf, n);
    csign_intarray_float((uint32_t *)zf, n);
    select_tests_det_float(yf, xf, zf, n, fpopts, MINMODE, MAXMODE,
                           i, i, 1, 0);
    free(xf);
    free(zf);
    free_array_float(yf, mode);
  }
}

/* Exactly representable normals. */
#test no_rounding_normal_numbers
printf("1c. No rounding: normal numbers\n");
for (size_t mode=1; mode<3; mode++) {
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    size_t n = ldexp(1., fpopts->precision - 1) * (ldexp(1., 5) - 2);
    uint64_t *xd = malloc(n * sizeof(*xd));
    init_intarray_double(xd, n, intminnormal_double(fpopts),
                         1ul << (52-fpopts->precision + 1));
    double *zd = alloc_init_array_double((double *)xd, n);
    double *yd = allocate_array_double((double *)xd, n, mode);
    select_tests_det_double(yd, (double *)xd, zd, n, fpopts,
                            MINMODE, MAXMODE, i, i, -1, -1);
    csign_intarray_double(xd, n);
    csign_intarray_double((uint64_t *)zd, n);
    select_tests_det_double(yd, (double *)xd, zd, n, fpopts,
                            MINMODE, MAXMODE, i, i, -1, -1);
    free(xd);
    free(zd);
    free_array_double(yd, mode);

    uint32_t *xf = malloc(n * sizeof(*xf));
    init_intarray_float(xf, n, intminnormal_float(fpopts),
                        1ul << (23-fpopts->precision + 1));
    float *zf = alloc_init_array_float((float *)xf, n);
    float *yf = allocate_array_float((float *)xf, n, mode);
    select_tests_det_float(yf, (float *)xf, zf, n, fpopts,
                           MINMODE, MAXMODE, i, i, -1, -1);
    csign_intarray_float(xf, n);
    csign_intarray_float((uint32_t *)zf, n);
    select_tests_det_float(yf, (float *)xf, (float *)xf, n, fpopts,
                           MINMODE, MAXMODE, i, i, -1, -1);
    free(xf);
    free(zf);
    free_array_float(yf, mode);
  }
}

/******************************************************************************/

/* 2. Deterministic rounding. */
#tcase deterministic_rounding

/* Underflow when subnormals are supported. */
#test deterministic_rounding_underflow_without_subnormals
printf("2a. Deterministic rounding: underflow without subnormals\n");
size_t n = 8;
double *xd = malloc(n * sizeof(*xd));
double *refzerod = calloc(n, sizeof(*refzerod));
double *refxmind = malloc(n * sizeof(*refxmind));
float *xf = malloc(n * sizeof(*xf));
float *refzerof = calloc(n, sizeof(*refzerof));
float *refxminf = malloc(n * sizeof(*refxminf));
for (size_t mode=1; mode<3; mode++) {
  double *yd = allocate_array_double(xd, n, mode);
  float *yf = allocate_array_float(xf, n, mode);
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    double xmin = minnormal(fpopts);
    double snxmin = minsubnormal(fpopts);
    double halfxmin = xmin / 2;
    double xd_imm [] = {nextafter(0, INFINITY),
                        nextafter(snxmin, 0),
                        snxmin,
                        nextafter(snxmin, INFINITY),
                        nextafter(halfxmin, 0),
                        halfxmin,
                        nextafter(halfxmin, INFINITY),
                        nextafter(xmin, 0)};
    double refd [] = {0, 0, 0, 0, 0, xmin, xmin, xmin};
    for (size_t j = 0; j < n; j++)
      refxmind[j] = xmin;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, -1, -1, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, 0, 1);
    refd[5] = 0;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 2, 2, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 3, 3, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 4, 4, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 7, 7, i, i, 0, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    csign_intarray_double((uint64_t *)refxmind, n);
    refd[5] = -xmin;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, -1, -1, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, 0, 1);
    refd[5] = 0;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 2, 2, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 3, 3, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 4, 4, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 7, 7, i, i, 0, 1);

    float xf_imm [] = {nextafterf(0, INFINITY),
                       nextafterf(snxmin, 0),
                       snxmin,
                       nextafterf(snxmin, INFINITY),
                       nextafterf(halfxmin, 0),
                       halfxmin,
                       nextafterf(halfxmin, INFINITY),
                       nextafterf(xmin, 0)};
    float reff [] = {0, 0, 0, 0, 0, xmin, xmin, xmin};
    for (size_t j = 0; j < n; j++)
      refxminf[j] = xmin;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, -1, -1, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, 0, 1);
    reff[5] = 0;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 2, 2, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 3, 3, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 4, 4, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 7, 7, i, i, 0, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    csign_intarray_float((uint32_t *)refxminf, n);
    reff[5] = -xmin;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, -1, -1, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, 0, 1);
    reff[5] = 0;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 2, 2, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 3, 3, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 4, 4, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 7, 7, i, i, 0, 1);
  }
  free_array_double(yd, mode);
  free_array_float(yf, mode);
}
free(xd);
free(refzerod);
free(refxmind);
free(xf);
free(refzerof);
free(refxminf);


/* Underflow when subnormals are not supported. */
#test deterministic_rounding_underflow_with_subnormals
printf("2b. Deterministic rounding: underflow with subnormals\n");
size_t n = 5;
double *xd = malloc(n * sizeof(*xd));
double *refzerod = calloc(n, sizeof(*refzerod));
double *refxmind = malloc(n * sizeof(*refxmind));
float *xf = malloc(n * sizeof(*xf));
float *refzerof = calloc(n, sizeof(*refzerof));
float *refxminf = malloc(n * sizeof(*refxminf));
for (size_t mode=1; mode<3; mode++) {
  double *yd = allocate_array_double(xd, n, mode);
  float *yf = allocate_array_float(xf, n, mode);
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    double xmin = minsubnormal(fpopts);
    double halfxmin = xmin / 2;
    double xd_imm [] = {nextafter(0, INFINITY),
                        nextafter(halfxmin, 0),
                        halfxmin,
                        nextafter(halfxmin, INFINITY),
                        nextafter(xmin, 0)};
    double refd [] = {0, 0, xmin, xmin, xmin};
    for (size_t j = 0; j < n; j++)
      refxmind[j] = xmin;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, -1, -1, i, i, 1, 1);
    refd[2] = 0;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 2, 2, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 3, 3, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 4, 4, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 7, 7, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    csign_intarray_double((uint64_t *)refxmind, n);
    refd[2] = -xmin;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, -1, -1, i, i, 1, 1);
    refd[2] = 0;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 2, 2, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 3, 3, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refzerod, n, fpopts, 4, 4, i, i, 1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmind, n, fpopts, 7, 7, i, i, 1, 1);

    float xf_imm [] = {nextafterf(0, INFINITY),
                       nextafterf(halfxmin, 0),
                       halfxmin,
                       nextafterf(halfxmin, INFINITY),
                       nextafterf(xmin, 0)};
    float reff [] = {0, 0, xmin, xmin, xmin};
    for (size_t j = 0; j < n; j++)
      refxminf[j] = xmin;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, -1, -1, i, i, 1, 1);
    reff[2] = 0;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 2, 2, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 3, 3, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 4, 4, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 7, 7, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    csign_intarray_float((uint32_t *)refxminf, n);
    reff[2] = -xmin;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, -1, -1, i, i, 1, 1);
    reff[2] = 0;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 2, 2, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 3, 3, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refzerof, n, fpopts, 4, 4, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxminf, n, fpopts, 7, 7, i, i, 1, 1);
  }
  free_array_double(yd, mode);
  free_array_float(yf, mode);
}
free(xd);
free(refzerod);
free(refxmind);
free(xf);
free(refzerof);
free(refxminf);

/*  Overflow. */
#test deterministic_rounding_overflow
printf("2c. Deterministic rounding: overflow\n");
size_t n = 4;
double *xd = malloc(n * sizeof(*xd));
double *refxmaxd = malloc(n * sizeof(*refxmaxd));
double *refinfd = malloc(n * sizeof(*refinfd));
float *xf = malloc(n * sizeof(*xf));
float *refxmaxf = malloc(n * sizeof(*refxmaxf));
float *refinff = malloc(n * sizeof(*refinff));
for (size_t mode=1; mode<3; mode++) {
  double *yd = allocate_array_double(xd, n, mode);
  float *yf = allocate_array_float(xf, n, mode);
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    double xmax = maxnormal(fpopts);
    double xbound = maxbound(fpopts);
    double xd_imm [] = {nextafter(xmax, INFINITY),
                        nextafter(xbound, 0),
                    xbound,
                    nextafter(xbound, INFINITY)};
    double refd [] = {xmax, xmax, inf_double(), inf_double()};
    for (size_t j = 0; j < n; j++) {
      refxmaxd[j] = xmax;
      refinfd[j] = inf_double();
    }
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, -1, -1, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, -1, 1);
    refd[2] = xmax;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, -1, 1);
    refd[2] = inf_double();
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refinfd, n, fpopts, 2, 2, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmaxd, n, fpopts, 3, 4, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmaxd, n, fpopts, 7, 7, i, i, -1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    csign_intarray_double((uint64_t *)refxmaxd, n);
    csign_intarray_double((uint64_t *)refinfd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, -1, -1, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, -1, 1);
    refd[2] = -xmax;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmaxd, n, fpopts, 2, 2, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refinfd, n, fpopts, 3, 3, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmaxd, n, fpopts, 4, 4, i, i, -1, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refxmaxd, n, fpopts, 7, 7, i, i, -1, 1);

    float xf_imm [] = {nextafterf(xmax, INFINITY),
                       nextafterf(xbound, 0),
                   xbound,
                   nextafterf(xbound, INFINITY)};
    float reff [] = {xmax, xmax, inf_float(), inf_float()};
    for (size_t j = 0; j < n; j++) {
      refxmaxf[j] = xmax;
      refinff[j] = inf_float();
    }
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, -1, -1, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, -1, 1);
    reff[2] = xmax;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, -1, 1);
    reff[2] = inf_float();
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refinff, n, fpopts, 2, 2, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxmaxf, n, fpopts, 3, 4, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxmaxf, n, fpopts, 7, 7, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    csign_intarray_float((uint32_t *)refxmaxf, n);
    csign_intarray_float((uint32_t *)refinff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, -1, -1, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, -1, 1);
    reff[2] = -xmax;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxmaxf, n, fpopts, 2, 2, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refinff, n, fpopts, 3, 3, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxmaxf, n, fpopts, 4, 4, i, i, -1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, refxmaxf, n, fpopts, 7, 7, i, i, -1, 1);
  }
  free_array_double(yd, mode);
  free_array_float(yf, mode);
}
free(xd);
free(refxmaxd);
free(refinfd);
free(xf);
free(refxmaxf);
free(refinff);

/* Rounding of numbers in the subnormal range. */
#test deterministic_rounding_subnormal_numbers
printf("2d. Deterministic rounding: subnormal numbers\n");
for (size_t mode=1; mode<3; mode++) {
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    size_t n = 3 * ldexp(1.,fpopts->precision-1) - 1;
    double *xd_imm = malloc(n * sizeof(*xd_imm));
    double *xd = malloc(n * sizeof(*xd));
    double *yd = allocate_array_double(xd, n, mode);
    double *refd = malloc(n * sizeof(*refd));
    double stepd = minsubnormal(fpopts);
    refd[0] = stepd;
    refd[1] = 2 * stepd;
    refd[2] = refd[1];
    for (size_t j = 3; j < n; j++)
      refd[j] = refd[j-3] + stepd;
    init_fparray_rounding_double(xd_imm, n, minsubnormal(fpopts), minsubnormal(fpopts));
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            -1, -1, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            -1, -1, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    for (size_t j = 1; j < n; j+=3)
      refd[j] -= stepd;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            0, 0, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            0, 0, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    for (size_t j = 1; j < n; j+=6)
      refd[j] += stepd;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            1, 1, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            1, 1, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    refd[0] = 2*stepd;
    for (size_t j = 1; j < n; j++)
      refd[j] = refd[j-1] + ((j % 3 == 0) ? stepd : 0);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            2, 2, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            3, 3, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    for (size_t j = 0; j < n; j++)
      refd[j] -= stepd;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            3, 4, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            2, 2, i, i, 1, 1);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            4, 4, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    refd[0] = stepd;
    for (size_t j = 1; j < n; j++)
      refd[j] = refd[j-1] + (((j-3) % 6 == 0) ? 2*stepd : 0);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            7, 7, i, i, 1, 1);
    csign_intarray_double((uint64_t *)xd_imm, n);
    csign_intarray_double((uint64_t *)refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, xd_imm, n);
    select_tests_det_double(yd, (double *)xd, (double *)refd, n, fpopts,
                            7, 7, i, i, 1, 1);
    free(xd_imm);
    free(xd);
    free_array_double(yd, mode);
    free(refd);

    float *xf_imm = malloc(n * sizeof(*xf_imm));
    float *xf = malloc(n * sizeof(*xf));
    float *yf = allocate_array_float(xf, n, mode);
    float *reff = malloc(n * sizeof(*reff));
    float stepf = minsubnormal(fpopts);
    reff[0] = stepf;
    reff[1] = 2 * stepf;
    reff[2] = reff[1];
    for (size_t j = 3; j < n; j++)
      reff[j] = reff[j-3] + stepf;
    init_fparray_rounding_float(xf_imm, n, minsubnormal(fpopts),
                                minsubnormal(fpopts));
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           -1, -1, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           -1, -1, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    for (size_t j = 1; j < n; j+=3)
      reff[j] -= stepf;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           0, 0, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           0, 0, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    for (size_t j = 1; j < n; j+=6)
      reff[j] += stepf;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           1, 1, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           1, 1, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    reff[0] = 2*stepf;
    for (size_t j = 1; j < n; j++)
      reff[j] = reff[j-1] + ((j % 3 == 0) ? stepf : 0);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf_imm, (float *)reff, n, fpopts,
                           2, 2, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           3, 3, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    for (size_t j = 0; j < n; j++)
      reff[j] -= stepf;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           3, 4, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           2, 2, i, i, 1, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           4, 4, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    reff[0] = stepf;
    for (size_t j = 1; j < n; j++)
      reff[j] = reff[j-1] + (((j-3) % 6 == 0) ? 2*stepf : 0);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           7, 7, i, i, 1, 1);
    csign_intarray_float((uint32_t *)xf_imm, n);
    csign_intarray_float((uint32_t *)reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, xf_imm, n);
    select_tests_det_float(yf, (float *)xf, (float *)reff, n, fpopts,
                           7, 7, i, i, 1, 1);
    free(xf_imm);
    free(xf);
    free_array_float(yf, mode);
    free(reff);
  }
}

/* Rounding of numbers in the normal range. */
#test deterministic_rounding_normal_numbers
printf("2e. Deterministic rounding: normal numbers\n");
for (size_t mode=1; mode<3; mode++) {
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    size_t n = 3 * (ldexp(1., fpopts->precision - 1) * (ldexp(1., 5) - 2) - 1);
    uint64_t *xd_imm = malloc(n * sizeof(*xd_imm));
    uint64_t *refd = malloc(n * sizeof(*refd));
    double *xd = malloc(n * sizeof(*xd));
    double *yd = allocate_array_double(xd, n, mode);
    uint64_t stepd = 1ul << (52-fpopts->precision + 1);
    init_intarray_rounding_double(xd_imm, n, intminnormal_double(fpopts), stepd);
    refd[0] = intminnormal_double(fpopts);
    refd[1] = intminnormal_double(fpopts) + stepd;
    refd[2] = intminnormal_double(fpopts) + stepd;
    for (size_t j = 3; j < n; j++)
      refd[j] = refd[j-3] + stepd;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            -1, -1, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            -1, -1, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    for (size_t j = 1; j < n; j+=3)
      refd[j] -= stepd;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            0, 0, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            0, 0, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    for (size_t j = 4; j < n; j+=6)
      refd[j] += stepd;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            1, 1, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            1, 1, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    refd[0] = intminnormal_double(fpopts) + stepd;
    for (size_t j = 1; j < n; j++)
      refd[j] = refd[j-1] + ((j % 3 == 0) ? stepd : 0);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            2, 2, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            3, 3, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    for (size_t j = 0; j < n; j++)
      refd[j] -= stepd;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            3, 4, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            2, 2, i, i, -1, -1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            4, 4, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    refd[0] = intminnormal_double(fpopts) + stepd;
    for (size_t j = 1; j < n; j++)
      refd[j] = refd[j-1] + ((j % 6 == 0) ? 2*stepd : 0);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            7, 7, i, i, -1, -1);
    csign_intarray_double(xd_imm, n);
    csign_intarray_double(refd, n);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)xd_imm, n);
    select_tests_det_double(yd, xd, (double *)refd, n, fpopts,
                            7, 7, i, i, -1, -1);
    free(xd_imm);
    free(refd);
    free(xd);
    free_array_double(yd, mode);

    uint32_t *xf_imm = malloc(n * sizeof(*xf_imm));
    uint32_t *reff = malloc(n * sizeof(*reff));
    float *xf = malloc(n * sizeof(*xf));
    float *yf = allocate_array_float(xf, n, mode);
    uint32_t stepf = 1u << (23-fpopts->precision + 1);
    init_intarray_rounding_float(xf_imm, n, intminnormal_float(fpopts), stepf);
    reff[0] = intminnormal_float(fpopts);
    reff[1] = intminnormal_float(fpopts) + stepf;
    reff[2] = intminnormal_float(fpopts) + stepf;
    for (size_t j = 3; j < n; j++)
      reff[j] = reff[j-3] + stepf;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           -1, -1, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           -1, -1, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    for (size_t j = 1; j < n; j+=3)
      reff[j] -= stepf;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           0, 0, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           0, 0, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    for (size_t j = 4; j < n; j+=6)
      reff[j] += stepf;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           1, 1, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           1, 1, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    reff[0] = intminnormal_float(fpopts) + stepf;
    for (size_t j = 1; j < n; j++)
      reff[j] = reff[j-1] + ((j % 3 == 0) ? stepf : 0);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           2, 2, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           3, 3, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    for (size_t j = 0; j < n; j++)
      reff[j] -= stepf;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           3, 4, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           2, 2, i, i, -1, -1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           4, 4, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    reff[0] = intminnormal_float(fpopts) + stepf;
    for (size_t j = 1; j < n; j++)
      reff[j] = reff[j-1] + ((j % 6 == 0) ? 2*stepf : 0);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           7, 7, i, i, -1, -1);
    csign_intarray_float(xf_imm, n);
    csign_intarray_float(reff, n);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)xf_imm, n);
    select_tests_det_float(yf, xf, (float *)reff, n, fpopts,
                           7, 7, i, i, -1, -1);
    free(xf_imm);
    free(reff);
    free(xf);
    free_array_float(yf, mode);
  }
}

/******************************************************************************/

/* 3. Stochastic rounding. */
#tcase stochastic_rounding

/* Underflow. */
#test stochastic_rounding_underflow
printf("3a. Stochastic rounding: underflow\n");
size_t subnormal, n = 3;
double *tmpdin = malloc(NREPS * sizeof(*tmpdin));
double *tmpdout = malloc(NREPS * sizeof(*tmpdout));
float *tmpfin = malloc(NREPS * sizeof(*tmpfin));
float *tmpfout = malloc(NREPS * sizeof(*tmpfout));
double proundequi = 0.50;
for (size_t mode=1; mode<3; mode++) {
  for(subnormal = 0; subnormal <= 1; subnormal++) {
    for(size_t i = 0; i < nformats; i++) {
      fpopts->precision = precision[i];
      fpopts->emax = emax[i];
      double xmin;
      if (subnormal)
        xmin = minsubnormal(fpopts);
      else
        xmin = minnormal(fpopts);
      double xd [] = {xmin/4, xmin/2, xmin/4*3};
      double proundprop [] = {0.75, 0.50, 0.25};
      select_tests_stoc_double(tmpdin, tmpdout, xd, proundprop,
                               n, fpopts, 5, i, i, subnormal, 1, mode);
      select_tests_stoc_double(tmpdin, tmpdout, xd, &proundequi,
                               n, fpopts, 6, i, i, subnormal, 1, mode);
      csign_intarray_double((uint64_t *)xd, n);
      invprob_array_double(proundprop, n);
      select_tests_stoc_double(tmpdin, tmpdout, xd, proundprop,
                               n, fpopts, 5, i, i, subnormal, 1, mode);
      select_tests_stoc_double(tmpdin, tmpdout, xd, &proundequi,
                               n, fpopts, 6, i, i, subnormal, 1, mode);
      invprob_array_double(proundprop, n);
      float xf [] = {xmin / 4, xmin / 2, xmin / 4 * 3};
      select_tests_stoc_float(tmpfin, tmpfout, xf, proundprop,
                              n, fpopts, 5, i, i, subnormal, 1, mode);
      select_tests_stoc_float(tmpfin, tmpfout, xf, &proundequi,
                              n, fpopts, 6, i, i, subnormal, 1, mode);
      csign_intarray_float((uint32_t *)xf, n);
      invprob_array_double(proundprop, n);
      select_tests_stoc_float(tmpfin, tmpfout, xf, proundprop,
                              n, fpopts, 5, i, i, subnormal, 1, mode);
      select_tests_stoc_float(tmpfin, tmpfout, xf, &proundequi,
                              n, fpopts, 6, i, i, subnormal, 1, mode);
    }
  }
}
free(tmpdin);
free(tmpdout);
free(tmpfin);
free(tmpfout);

/* Overflow. */
#test stochastic_rounding_overflow
printf("3b. Stochastic rounding: overflow\n");
size_t subnormal, i, n = 3;
double *tmpdin = malloc(NREPS * sizeof(*tmpdin));
double *tmpdout = malloc(NREPS * sizeof(*tmpdout));
float *tmpfin = malloc(NREPS * sizeof(*tmpfin));
float *tmpfout = malloc(NREPS * sizeof(*tmpfout));
double proundequi = 0.50;
for (size_t mode=1; mode<3; mode++ ) {
  for(subnormal = 0; subnormal <= 1; subnormal++) {
    for(i = 0; i < nformats; i++) {
      fpopts->precision = precision[i];
      fpopts->emax = emax[i];
      double xmax = maxnormal(fpopts);
      double xbound = maxbound(fpopts);
      double *yd =  malloc(n * sizeof(*yd));
      double xdin [] = {xbound, nextafter(xbound, INFINITY), 2 * xbound};
      double xdout [] = {inf_double(), inf_double(), inf_double()};
      select_tests_det_double(yd, xdin, xdout, n, fpopts, 5, 5,
                              i, i, 1, 1);
      select_tests_stoc_double(tmpdin, tmpdout, xdin, &proundequi,
                               n, fpopts, 6, i, i, subnormal, 1, mode);
      csign_intarray_double((uint64_t *)xdin, n);
      csign_intarray_double((uint64_t *)xdout, n);
      select_tests_det_double(yd, xdin, xdout, n, fpopts, 5, 5,
                              i, i, 1, 1);
      select_tests_stoc_double(tmpdin, tmpdout, xdin, &proundequi,
                               n, fpopts, 6, i, i, subnormal, 1, mode);
      double xd [] = {(3*xmax+xbound) / 4,
                      (xmax+xbound) / 2,
                      (xmax+3*xbound) / 4};
      double proundprop [] = {0.75, 0.5, 0.25};
      select_tests_stoc_double(tmpdin, tmpdout, xd, proundprop,
                               n, fpopts, 5, i, i, subnormal, 1, mode);
      select_tests_stoc_double(tmpdin, tmpdout, xd, &proundequi,
                               n, fpopts, 6, i, i, subnormal, 1, mode);
      csign_intarray_double((uint64_t *)xd, n);
      invprob_array_double(proundprop, n);
      select_tests_stoc_double(tmpdin, tmpdout, xd, proundprop,
                               n, fpopts, 5, i, i, subnormal, 1, mode);
      select_tests_stoc_double(tmpdin, tmpdout, xd, &proundequi,
                               n, fpopts, 6, i, i, subnormal, 1, mode);
      invprob_array_double(proundprop, n);
      float *yf =  malloc(n * sizeof(*yf));
      if (i == 0) {
        float xfin [] = {xbound, nextafterf(xbound, INFINITY), 2 * xbound};
        float xfout [] = {inf_double(), inf_double(), inf_double()};
        select_tests_det_float(yf, xfin, xfout, n, fpopts, 5, 5,
                               i, i, 1, 1);
        select_tests_stoc_float(tmpfin, tmpfout, xfin, &proundequi,
                                n, fpopts, 6, i, i, subnormal, 1, mode);
        csign_intarray_float((uint32_t *)xfin, n);
        csign_intarray_float((uint32_t *)xfout, n);
        select_tests_det_float(yf, xfin, xfout, n, fpopts, 5, 5,
                               i, i, 1, 1);
        select_tests_stoc_float(tmpfin, tmpfout, xfin, &proundequi,
                                n, fpopts, 6, i, i, subnormal, 1, mode);
        float xf [] = {(3*xmax+xbound) / 4,
                       (xmax+xbound) / 2,
                       (xmax+3*xbound) / 4};
        select_tests_stoc_float(tmpfin, tmpfout, xf, proundprop,
                                n, fpopts, 5, i, i, subnormal, 1, mode);
        select_tests_stoc_float(tmpfin, tmpfout, xf, &proundequi,
                                n, fpopts, 6, i, i, subnormal, 1, mode);
        csign_intarray_float((uint32_t *)xf, n);
        invprob_array_double(proundprop, n);
        select_tests_stoc_float(tmpfin, tmpfout, xf, proundprop,
                                n, fpopts, 5, i, i, subnormal, 1, mode);
        select_tests_stoc_float(tmpfin, tmpfout, xf, &proundequi,
                                n, fpopts, 6, i, i, subnormal, 1, mode);
        invprob_array_double(proundprop, n);
      }
    }
  }
}
free(tmpdin);
free(tmpdout);
free(tmpfin);
free(tmpfout);

/* Rounding of numbers in the subnormal range. */
#test stochastic_rounding_subnormal_numbers
printf("3c. Stochastic rounding: subnormal numbers\n");
double *tmpdin = malloc(NREPS * sizeof(*tmpdin));
double *tmpdout = malloc(NREPS * sizeof(*tmpdout));
float *tmpfin = malloc(NREPS * sizeof(*tmpfin));
float *tmpfout = malloc(NREPS * sizeof(*tmpfout));
double proundequi = 0.50;
for (size_t mode=1; mode<3; mode++ ) {
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];
    size_t n = 3 * ldexp(1.,fpopts->precision-1) - 1;
    double *xd = malloc(n * sizeof(*xd));
    double stepd = minsubnormal(fpopts);
    init_fparray_rounding_stoc_double(xd, n, stepd, stepd);
    double proundprop [] = {0.75, 0.50, 0.25};
    select_tests_stoc_double(tmpdin, tmpdout, xd, proundprop,
                             n, fpopts, 5, i, i, 1, 1, mode);
    select_tests_stoc_double(tmpdin, tmpdout, xd, &proundequi,
                             n, fpopts, 6, i, i, 1, 1, mode);
    csign_intarray_double((uint64_t *)xd, n);
    invprob_array_double(proundprop, 3);
    select_tests_stoc_double(tmpdin, tmpdout, xd, proundprop,
                             n, fpopts, 5, i, i, 1, 1, mode);
    select_tests_stoc_double(tmpdin, tmpdout, xd, &proundequi,
                             n, fpopts, 6, i, i, 1, 1, mode);
    invprob_array_double(proundprop, 3);
    free(xd);

    float *xf = malloc(n * sizeof(*xf));
    float stepf = minsubnormal(fpopts);
    init_fparray_rounding_stoc_float(xf, n, stepf, stepf);
    select_tests_stoc_float(tmpfin, tmpfout, xf, proundprop,
                            n, fpopts, 5, i, i, 1, 1, mode);
    select_tests_stoc_float(tmpfin, tmpfout, xf, &proundequi,
                            n, fpopts, 6, i, i, 1, 1, mode);
    csign_intarray_float((uint32_t *)xf, n);
    invprob_array_double(proundprop, 3);
    select_tests_stoc_float(tmpfin, tmpfout, xf, proundprop,
                            n, fpopts, 5, i, i, 1, 1, mode);
    select_tests_stoc_float(tmpfin, tmpfout, xf, &proundequi,
                            n, fpopts, 6, i, i, 1, 1, mode);
    free(xf);
  }
 }
free(tmpdin);
free(tmpdout);
free(tmpfin);
free(tmpfout);

/* Rounding of numbers in the normal range. */
#test stochastic_rounding_normal_numbers
printf("3d. Stochastic rounding: normal numbers\n");
size_t subnormal, explim;
double *tmpdin = malloc(NREPS * sizeof(*tmpdin));
double *tmpdout = malloc(NREPS * sizeof(*tmpdout));
float *tmpfin = malloc(NREPS * sizeof(*tmpfin));
float *tmpfout = malloc(NREPS * sizeof(*tmpfout));
double proundequi = 0.50;
for (size_t mode=1; mode<3; mode++ ) {
  for (subnormal = 0; subnormal <= 1; subnormal++) {
    for (explim = 0; explim <= 1; explim++) {
      for (size_t i = 0; i < nformats; i++) {
        fpopts->precision = precision[i];
        fpopts->emax = emax[i];
        size_t n = 3 * (ldexp(1., fpopts->precision - 1) * (ldexp(1., 5) - 2) - 1);
        uint64_t *xd = malloc(n * sizeof(*xd));
        double xmin = minnormal(fpopts);
        uint64_t stepd = 1ul << (52-fpopts->precision + 1);
        init_intarray_rounding_stoc_double(xd, n, INTOFd(xmin), stepd);
        double proundprop [] = {0.75, 0.50, 0.25};
        select_tests_stoc_double(tmpdin, tmpdout, (double *)xd, proundprop,
                                 n, fpopts, 5, i, i, subnormal, explim, mode);
        select_tests_stoc_double(tmpdin, tmpdout, (double *)xd, &proundequi,
                                 n, fpopts, 6, i, i, subnormal, explim, mode);
        csign_intarray_double(xd, n);
        invprob_array_double(proundprop, 3);
        select_tests_stoc_double(tmpdin, tmpdout, (double *)xd, proundprop,
                                 n, fpopts, 5, i, i, subnormal, explim, mode);
        select_tests_stoc_double(tmpdin, tmpdout, (double *)xd, &proundequi,
                                 n, fpopts, 6, i, i, subnormal, explim, mode);
        invprob_array_double(proundprop, 3);
        free(xd);
        uint32_t *xf = malloc(n * sizeof(*xf));
        uint32_t stepf = 1ul << (23-fpopts->precision + 1);
        float xminf = minnormal(fpopts);
        init_intarray_rounding_stoc_float(xf, n, INTOFf(xminf), stepf);
        select_tests_stoc_float(tmpfin, tmpfout, (float *)xf, proundprop,
                                n, fpopts, 5, i, i, subnormal, explim, mode);
        select_tests_stoc_float(tmpfin, tmpfout, (float *)xf, &proundequi,
                                n, fpopts, 6, i, i, subnormal, explim, mode);
        csign_intarray_float((uint32_t *)xf, n);
        invprob_array_double(proundprop, 3);
        select_tests_stoc_float(tmpfin, tmpfout, (float *)xf, proundprop,
                                n, fpopts, 5, i, i, subnormal, explim, mode);
        select_tests_stoc_float(tmpfin, tmpfout, (float *)xf, &proundequi,
                                n, fpopts, 6, i, i, subnormal, explim, mode);
        free(xf);
      }
    }
  }
}
free(tmpdin);
free(tmpdout);
free(tmpfin);
free(tmpfout);

#main-pre
putenv("CK_VERBOSITY=verbose");
tcase_add_unchecked_fixture(tc1_1, fpopts_setup, fpopts_teardown);
tcase_add_unchecked_fixture(tc1_2, fpopts_setup, fpopts_teardown);
tcase_add_unchecked_fixture(tc1_3, fpopts_setup, fpopts_teardown);
tcase_set_timeout(tc1_1, 100);
tcase_set_timeout(tc1_2, 200);
tcase_set_timeout(tc1_3, 300);
srunner_set_log(sr, "cpfloat_test.log");

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
