/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

#include <check.h>
#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <inttypes.h>
#include "cpfloat_binary32.h"
#include "cpfloat_binary64.h"

/* Struct to reinterpret floating-point values as unsigned integers. */
typedef union {
  uint32_t intval;
  float fpval;
} fpint_f;

#define INTOFf(x)(((fpint_f){.fpval = (float)(x)}).intval)
#define FPOFf(x)(((fpint_f){.intval = (uint32_t)(x)}).fpval)

typedef union {
  uint64_t intval;
  double fpval;
} fpint_d;

#define INTOFd(x)(((fpint_d){.fpval = (double)(x)}).intval)
#define FPOFd(x)(((fpint_d){.intval = (uint64_t)(x)}).fpval)



#define CONST_E         2.7182818284590452354
#define CONST_LOG2E     1.4426950408889634074
#define CONST_PI        3.14159265358979323846
#define CONST_SQRT2     1.41421356237309504880

#define MINFORMAT 0
#define MAXFORMAT 2
#define MINMODE  -1
#define MAXMODE   8
#define NREPS  1000

/* Define default target formats. */
// binary16, bfloat16, TensorFloat-32
static size_t precision [] = {11, 8, 11};
static size_t emax [] = {15, 127, 127};
static size_t nformats = 2;

/* Structure for options and fixtures. */
optstruct *fpopts;

void fpopts_setup(void) {
  fpopts = malloc(sizeof(optstruct));
  fpopts->flip = 0;
  fpopts->p = 0.5;
  fpopts->round = 1;
  fpopts->explim = CPFLOAT_EXPRANGE_STOR;
  fpopts->bitseed = NULL;
  fpopts->randseedf = NULL;
  fpopts->randseed = NULL;
}

void fpopts_teardown(void) {
  free(fpopts);
}

/***********************
 * Auxiliary functions *
 ***********************/

/* Return values of interest */
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
uint64_t intmaxsubnormal_double(optstruct *fpopts) {
  return INTOFd(ldexp(1., 1 - fpopts->emax) *
               (1 - ldexp(1., 1 - fpopts->precision)));
}

static inline
uint32_t intmaxsubnormal_float(optstruct *fpopts) {
  return INTOFf(ldexp(1., 1 - fpopts->emax) *
               (1 - ldexp(1., 1 - fpopts->precision)));
}

static inline
double minnormal(optstruct *fpopts) {
  return ldexp(1., 1-fpopts->emax);
}

static inline
uint64_t intminnormal_double(optstruct *fpopts) {
  return INTOFd(ldexp(1., 1-fpopts->emax));
}

static inline
uint32_t intminnormal_float(optstruct *fpopts) {
  return INTOFf(ldexp(1., 1-fpopts->emax));
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

/* Initialize arrays. */
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

/* Change sign of arrays. */
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


/* Compare result of Booelan functions. */
static inline
bool both_true_or_false(bool a, bool b) {
  return !(a ^ b);
}

static inline
void check_equality_double(double *x, double *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_double(x[j], y[j])) {
      printf("DOUBLE\n");
      printf("***\nj = %ld\nx  = %23.15e [%" PRIu64 "]\ny = %23.15e [%" PRIu64 "]\n",
             j,
             x[j], *(uint64_t *)(x + j),
             y[j], *(uint64_t *)(y + j));
    }
    ck_assert(nan_safe_compare_double(x[j], y[j]));
  }
}

static inline
void check_equality_double_int(double *x, int *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_double(x[j], y[j])) {
      printf("DOUBLE\n");
      printf("***\nj = %ld\nx  = %23.15e [%" PRIu64 "]\ny = %23.15e\n",
             j, x[j], *(uint64_t *)(x + j), (double)y[j]);
    }
    ck_assert(nan_safe_compare_double(x[j], (double)y[j]));
  }
}

static inline
void check_equality_double_long(double *x, long *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_double(x[j], y[j])) {
      printf("DOUBLE\n");
      printf("***\nj = %ld\nx  = %23.15e [%" PRIu64 "]\ny = %23.15e\n",
             j, x[j], *(uint64_t *)(x + j), (double)y[j]);
    }
    ck_assert(nan_safe_compare_double(x[j], (double)y[j]));
  }
}

static inline
void check_equality_double_long_long(double *x, long long *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_double(x[j], y[j])) {
      printf("DOUBLE\n");
      printf("***\nj = %ld\nx  = %23.15e [%" PRIu64 "]\ny = %23.15e\n",
             j, x[j], *(uint64_t *)(x + j), (double)y[j]);
    }
    ck_assert(nan_safe_compare_double(x[j], (double)y[j]));
  }
}

static inline
void check_equality_float(float *x, float *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
        if (!nan_safe_compare_float(x[j], y[j])) {
      printf("FLOAT\n");
      printf("***\nj = %ld\nx  = %23.15e [%X]\ny = %23.15e\n",
             j, x[j], *(uint32_t *)(x + j), (float)y[j]);
    }
    ck_assert(nan_safe_compare_float(x[j], y[j]));
  }
}

static inline
void check_equality_float_int(float *x, int *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_float(x[j], y[j])) {
      printf("FLOAT\n");
      printf("***\nj = %ld\nx  = %23.15e [%X]\ny = %23.15e\n",
             j, x[j], *(uint32_t *)(x + j), (float)y[j]);
    }
    ck_assert(nan_safe_compare_float(x[j], (float)y[j]));
  }
}

static inline
void check_equality_float_long(float *x, long *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_float(x[j], y[j])) {
      printf("FLOAT\n");
      printf("***\nj = %ld\nx  = %23.15e [%X]\ny = %23.15e\n",
             j, x[j], *(uint32_t *)(x + j), (float)y[j]);
    }
    ck_assert(nan_safe_compare_float(x[j], (float)y[j]));
  }
}

static inline
void check_equality_float_long_long(float *x, long long *y, size_t n) {
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_float(x[j], y[j])) {
      printf("FLOAT\n");
      printf("***\nj = %ld\nx  = %23.15e [%X]\ny = %23.15e\n",
             j, x[j], *(uint32_t *)(x + j), (float)y[j]);
    }
    ck_assert(nan_safe_compare_float(x[j], (float)y[j]));
  }
}


/* Run on single array and check results. */
static inline
void check_array_double(double *y, double *x, double *ref,
                        size_t n, optstruct * fpopts) {
  cpfloat(y, x, n, fpopts);
  for (size_t j = 0; j < n; j++) {
    if (!nan_safe_compare_double(y[j], ref[j])) {
      printf("DOUBLE\n");
      printf("***\n"
             "j = %ld\n"
             "in  = %23.15e [%" PRIu64 "]\n"
             "out = %23.15e [%" PRIu64 "]\n"
             "ref = %23.15e [%" PRIu64 "]\n"
             "r = %d\n"
             "pr = %d\n"
             "e = %d\n"
             "s = %d\n",
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
      printf("***\n"
             "j = %ld\n"
             "in  = %16.8e [%X]\n"
             "out = %16.8e [%X]\n"
             "ref = %16.8e [%X]\n"
             "r = %d\n"
             "pr = %d\n"
             "e = %d\n"
             "s = %d\n",
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
      printf("%e\n", fabs(counter[0]/NREPS - prounddown[i % 3]));
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
      printf("%e\n", fabs(counter[0]/NREPS - prounddown[i % 3]));
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
      printf("%e\n", fabs(counter[0]/NREPS - prounddown[i % 3]));
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
      printf("%e\n", fabs(counter[0]/NREPS - prounddown[i % 3]));
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

double ref_add(double a, double b) {return a + b;}
double ref_sub(double a, double b) {return a - b;}
double ref_mul(double a, double b) {return a * b;}
double ref_div(double a, double b) {return a / b;}

float ref_addf(float a, float b) {return a + b;}
float ref_subf(float a, float b) {return a - b;}
float ref_mulf(float a, float b) {return a * b;}
float ref_divf(float a, float b) {return a / b;}

void test_univariate_arith_op_double(double* xd,
                              double *ref,
                              const double *a,
                              const size_t n,
                              optstruct *fpopts,
                              int (*cpf_fun)(double *,
                                             const double *,
                                             const size_t,
                                             optstruct *),
                              double (*fun)(double)) {
  for (size_t i = 0; i < n; i++)
    ref[i] = fun(a[i]);
  ck_assert(cpf_fun(xd, a, n, fpopts) == 0);
  check_array_double(ref, ref, xd, n, fpopts);
}

void test_bivariate_arith_op_double(double* xd,
                             double *ref,
                             const double *a,
                             const double *b,
                             const size_t n,
                             optstruct *fpopts,
                             int (*cpf_fun)(double *,
                                            const double *,
                                            const double *,
                                            const size_t,
                                            optstruct *),
                             double (*fun)(double, double)) {
  for (size_t i = 0; i < n; i++)
    ref[i] = fun(a[i], b[i]);
  ck_assert(cpf_fun(xd, a, b, n, fpopts) == 0);
  check_array_double(ref, ref, xd, n, fpopts);
}

void test_trivariate_arith_op_double(double* xd,
                              double *ref,
                              const double *a,
                              const double *b,
                              const double *c,
                              const size_t n,
                              optstruct *fpopts,
                              int (*cpf_fun)(double *,
                                             const double *,
                                             const double *,
                                             const double *,
                                             const size_t,
                                             optstruct *),
                              double (*fun)(double, double, double)) {
  for (size_t i = 0; i < n; i++)
    ref[i] = fun(a[i], b[i], c[i]);
  ck_assert(cpf_fun(xd, a, b, c, n, fpopts) == 0);
  check_array_double(ref, ref, xd, n, fpopts);
}

void test_univariate_arith_op_float(float* xd,
                                    float *ref,
                                    const float *a,
                                    const size_t n,
                                    optstruct *fpopts,
                                    int (*cpf_fun)(float *,
                                                   const float *,
                                                   const size_t,
                                                   optstruct *),
                                    float (*fun)(float)) {
  for (size_t i = 0; i < n; i++)
    ref[i] = fun(a[i]);
  ck_assert(cpf_fun(xd, a, n, fpopts) == 0);
  check_array_float(ref, ref, xd, n, fpopts);
}

void test_bivariate_arith_op_float(float* xd,
                                   float *ref,
                                   const float *a,
                                   const float *b,
                                   const size_t n,
                                   optstruct *fpopts,
                                   int (*cpf_fun)(float *,
                                                  const float *,
                                                  const float *,
                                                  const size_t,
                                                  optstruct *),
                                   float (*fun)(float, float)) {
  for (size_t i = 0; i < n; i++)
    ref[i] = fun(a[i], b[i]);
  ck_assert(cpf_fun(xd, a, b, n, fpopts) == 0);
  check_array_float(ref, ref, xd, n, fpopts);
}

void test_trivariate_arith_op_float(float* xd,
                                    float *ref,
                                    const float *a,
                                    const float *b,
                                    const float *c,
                                    const size_t n,
                                    optstruct *fpopts,
                                    int (*cpf_fun)(float *,
                                                   const float *,
                                                   const float *,
                                                   const float *,
                                                   const size_t,
                                                   optstruct *),
                                    float (*fun)(float, float, float)) {
  for (size_t i = 0; i < n; i++)
    ref[i] = fun(a[i], b[i], c[i]);
  ck_assert(cpf_fun(xd, a, b, c, n, fpopts) == 0);
  check_array_float(ref, ref, xd, n, fpopts);
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
    size_t n = ldexp(1., fpopts->precision-1) * 2 * fpopts->emax;
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
    refd[5] = 0;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, 0, 1);
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
    refd[5] = 0;
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_double(yd, n, mode);
    copy_array_double(xd, (double *)&xd_imm, n);
    select_tests_det_double(yd, xd, refd, n, fpopts, 1, 1, i, i, 0, 1);
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
    reff[5] = 0;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, 0, 1);
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
    reff[5] = 0;
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 0, 0, i, i, 0, 1);
    reset_array_float(yf, n, mode);
    copy_array_float(xf, (float *)&xf_imm, n);
    select_tests_det_float(yf, xf, reff, n, fpopts, 1, 1, i, i, 0, 1);
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
    init_fparray_rounding_double(xd_imm, n,
                                 minsubnormal(fpopts), minsubnormal(fpopts));
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
    size_t n = 3 * (ldexp(1., fpopts->precision-1) * 2 * fpopts->emax - 1);
    uint64_t *xd_imm = malloc(n * sizeof(*xd_imm));
    uint64_t *refd = malloc(n * sizeof(*refd));
    double *xd = malloc(n * sizeof(*xd));
    double *yd = allocate_array_double(xd, n, mode);
    uint64_t stepd = 1ul << (52-fpopts->precision + 1);
    init_intarray_rounding_double(xd_imm, n,
                                  intminnormal_double(fpopts), stepd);
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
    init_intarray_rounding_float(xf_imm, n,
                                 intminnormal_float(fpopts), stepf);
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

        size_t n = 3 * (ldexp(1., fpopts->precision-1) * 2 * fpopts->emax - 1);
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

/******************************************************************************/

/* 4. Floating-point operations. */
#tcase operations

#test arithmetic_operations_correctness
printf("4a. Correctness of arithmetic operations\n");
double xmin = minnormal(fpopts);
double xmins = minsubnormal(fpopts);
double ad[] = {nan_double(), -inf_double(), -2., -1, -xmin, -xmins, 0.,
               xmins, xmin, 1., CONST_SQRT2, CONST_LOG2E, 2., CONST_E,
               CONST_PI, inf_double()};
double bd[] = {nan_double(), -inf_double(), -2., -1, -xmin, -xmins, 0.,
               xmins, xmin, 1., CONST_SQRT2, CONST_LOG2E, 2., CONST_E,
               CONST_PI, inf_double()};
double cd[] = {nan_double(), -inf_double(), -2., -1, -xmin, -xmins, 0.,
               xmins, xmin, 1., CONST_SQRT2, CONST_LOG2E, 2., CONST_E,
               CONST_PI, inf_double()};

size_t n = 16;
double *xd = malloc(n * sizeof(*xd));
double *refd = malloc(n * sizeof(* xd));

int det_roundings[] = {-1, 0, 1, 2, 3, 4, 7, 8};
size_t n_roundings = 7;

for (size_t i = 0; i < nformats; i++) {
  fpopts->precision = precision[i];
  fpopts->emax = emax[i];
  for (size_t i = 0; i < n_roundings; i++) {
    fpopts->round = det_roundings[i];

    // Univariate functions.
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_cos, cos);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_sin, sin);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_tan, tan);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_acos, acos);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_asin, asin);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_atan, atan);

    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_cosh, cosh);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_sinh, sinh);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_tanh, tanh);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_acosh, acosh);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_asinh, asinh);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_atanh, atanh);

    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_exp, exp);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_log, log);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_log10, log10);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_exp2, exp2);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_expm1, expm1);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_log1p, log1p);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_log2, log2);

    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_sqrt, sqrt);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_cbrt, cbrt);

    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_erf, erf);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_erfc, erfc);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_tgamma, tgamma);
    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_lgamma, lgamma);

    double *ones = malloc(n * sizeof(* ones));
    for (size_t j = 0; j < n; j++) {
      ones[j] = 1.;
      refd[j] = ad[j] < 0 ? -1. : 1.;
    }
    test_bivariate_arith_op_double(xd, refd, ones, ad, n, fpopts,
                                   cpf_copysign, copysign);
    free(ones);

    test_univariate_arith_op_double(xd, refd, ad, n, fpopts,
                                    cpf_fabs, fabs);

    for (size_t j =0; j < n; j++) {
      // Bivariate functions.
      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_add, ref_add);
      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_sub, ref_sub);
      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_mul, ref_mul);
      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_div, ref_div);

      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_atan2, atan2);

      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_pow, pow);

      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_hypot, hypot);

      test_bivariate_arith_op_double(xd, refd, ad, bd, n, fpopts,
                                     cpf_fmod, fmod);

      // Trivariate functions.
      for(size_t k = 0; k < n; k++) {
        test_trivariate_arith_op_double(xd, refd,
                                        ad, bd, cd, n, fpopts,
                                        cpf_fma, fma);
        double cfirst  = cd[0];
        for (size_t l = 0; l < n-1; l++)
          cd[l] = cd[l+1];
        cd[n-1] = cfirst;
      }

      double bfirst = bd[0];
      for (size_t l = 0; l < n-1; l++)
        bd[l] = bd[l+1];
      bd[n-1] = bfirst;
    }
  }
}
free(xd);
free(refd);

float xminf = minnormal(fpopts);
float xminfs = minsubnormal(fpopts);
float af[] = {nan_float(), -inf_float(), -2., -1, -xminf, -xminfs, 0.,
               xminfs, xminf, 1., CONST_SQRT2, CONST_LOG2E, 2., CONST_E,
               CONST_PI, inf_float()};
float bf[] = {nan_float(), -inf_float(), -2., -1, -xminf, -xminfs, 0.,
               xminfs, xminf, 1., CONST_SQRT2, CONST_LOG2E, 2., CONST_E,
               CONST_PI, inf_float()};
float cf[] = {nan_float(), -inf_float(), -2., -1, -xminf, -xminfs, 0.,
               xminfs, xminf, 1., CONST_SQRT2, CONST_LOG2E, 2., CONST_E,
               CONST_PI, inf_float()};

float *xf = malloc(n * sizeof(*xf));
float *reff = malloc(n * sizeof(* xf));

for (size_t i = 0; i < nformats; i++) {
  fpopts->precision = precision[i];
  fpopts->emax = emax[i];
  for (size_t i = 0; i < n_roundings; i++) {
    fpopts->round = det_roundings[i];
    // Univariate functions.
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_cosf, cosf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_sinf, sinf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_tanf, tanf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_acosf, acosf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_asinf, asinf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_atanf, atanf);

    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_coshf, coshf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_sinhf, sinhf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_tanhf, tanhf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_acoshf, acoshf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_asinhf, asinhf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_atanhf, atanhf);

    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_expf, expf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_logf, logf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_log10f, log10f);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_exp2f, exp2f);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_expm1f, expm1f);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_log1pf, log1pf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_log2f, log2f);

    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_sqrtf, sqrtf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_cbrtf, cbrtf);

    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_erff, erff);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_erfcf, erfcf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_tgammaf, tgammaf);
    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_lgammaf, lgammaf);

    float *ones = malloc(n * sizeof(* ones));
    for (size_t j = 0; j < n; j++) {
      ones[j] = 1.;
      reff[j] = af[j] < 0 ? -1. : 1.;
    }
    test_bivariate_arith_op_float(xf, reff, ones, af, n, fpopts,
                                   cpf_copysignf, copysignf);
    free(ones);

    test_univariate_arith_op_float(xf, reff, af, n, fpopts,
                                   cpf_fabsf, fabsf);

    for (size_t j =0; j < n; j++) {
      // Bivariate functions.
      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_addf, ref_addf);
      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_subf, ref_subf);
      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_mulf, ref_mulf);
      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_divf, ref_divf);

      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_atan2f, atan2f);

      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_powf, powf);

      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_hypotf, hypotf);

      test_bivariate_arith_op_float(xf, reff, af, bf, n, fpopts,
                                    cpf_fmodf, fmodf);

      // Trivariate functions.
      for(size_t k = 0; k < n; k++) {
        test_trivariate_arith_op_float(xf, reff,
                                       af, bf, cf, n, fpopts,
                                       cpf_fmaf, fmaf);
        float cfirst  = cf[0];
        for (size_t l = 0; l < n-1; l++)
          cf[l] = cf[l+1];
        cf[n-1] = cfirst;
      }
      float bfirst = bf[0];
      for (size_t l = 0; l < n-1; l++)
        bf[l] = bf[l+1];
      bf[n-1] = bfirst;
    }
  }
}
free(xf);
free(reff);

#test floating_point_operations
printf("4b. Floating-point operations\n");
size_t n = 16;
int *exp = malloc(n * sizeof(*exp));
long *lexp = malloc(n * sizeof(*lexp));

double *xd = malloc(n * sizeof(*xd));
double *intpartd = malloc(n * sizeof(*intpartd));
double *refd = malloc(n * sizeof(*refd));
double *oned = malloc(n * sizeof(*oned));
for (size_t j = 0; j < n; j++)
  oned[j] = 1.0;

int det_roundings[] = {-1, 0, 1, 2, 3, 4, 7, 8};
size_t n_roundings = 7;

const int resd_isfinite [] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0};
const int resd_isinf    [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};
const int resd_isnan    [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
const int resd_isnormal [] = {1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
const int resd_fpclassify [] = {FP_NORMAL, FP_NORMAL,FP_NORMAL, FP_SUBNORMAL,
                                FP_SUBNORMAL, FP_NORMAL, FP_NORMAL, FP_NORMAL,
                                FP_NORMAL, FP_NORMAL, FP_NORMAL, FP_NORMAL,
                                FP_NAN, FP_ZERO, FP_INFINITE, FP_INFINITE};

for (size_t i = 0; i < nformats; i++) {
  fpopts->precision = precision[i];
  fpopts->emax = emax[i];
  double xmind = minnormal(fpopts);
  double xminds = minsubnormal(fpopts);
  double ad[] = {-2., -1, -xmind, -xminds, xminds, xmind, 1., CONST_SQRT2,
                 CONST_LOG2E, 2., CONST_E, CONST_PI,
                 nan_double(), 0, -inf_double(), inf_double()};

  for (size_t i = 0; i < n_roundings; i++) {
    fpopts->round = det_roundings[i];

    cpfloat(refd, ad, n, fpopts);
    cpf_frexp(xd, exp, ad, n, fpopts);
    cpf_ldexp(xd, xd, exp, n, fpopts);
    check_equality_double(xd, refd, n);

    cpf_modf(xd, intpartd, ad, n, fpopts);
    for (size_t j = 0; j < n; j++)
      xd[j] = xd[j] + intpartd[j];
    check_equality_double(xd, refd, n);

    // These test would not work for +inf, -inf, and 0 because of the division.
    cpf_ilogb(exp, ad, n-3, fpopts);
    cpf_scalbn(xd, oned, exp, n-3, fpopts);
    cpf_div(xd, ad, xd, n-3, fpopts);
    cpf_ldexp(xd, xd, exp, n-3, fpopts);
    check_equality_double(xd, refd, n-3);

    cpf_ilogb(exp, ad, n-3, fpopts);
    for (size_t j = 0; j < n; j++)
      lexp[j] = (long)exp[j];
    cpf_scalbln(xd, oned, lexp, n-3, fpopts);
    cpf_div(xd, ad, xd, n-3, fpopts);
    cpf_ldexp(xd, xd, exp, n-3, fpopts);
    check_equality_double(xd, refd, n-3);

    // Removing nans and infs for this test.
    cpf_ilogb(exp, ad, n-4, fpopts);
    cpf_logb(xd, ad, n-4, fpopts);
    for (size_t j = 0; j < n-4; j++)
      ck_assert(nan_safe_compare_double(xd[j], (double) exp[j]));

    cpf_isfinite(exp, ad, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(both_true_or_false(exp[j], resd_isfinite[j]));

    cpf_isinf(exp, ad, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(both_true_or_false(exp[j], resd_isinf[j]));

    cpf_isnan(exp, ad, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(both_true_or_false(exp[j], resd_isnan[j]));

    // Subnormals depend on the target format: we have to use the correct
    // exponent range and the correct subnormal handling mode to get the
    // reference results.
    fpopts->explim = CPFLOAT_EXPRANGE_TARG;
    fpopts->subnormal = CPFLOAT_SUBN_USE;
    cpf_isnormal(exp, ad, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(both_true_or_false(exp[j], resd_isnormal[j]));
    cpf_fpclassify(exp, ad, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(exp[j] == resd_fpclassify[j]);
    fpopts->subnormal = CPFLOAT_SUBN_RND;
    fpopts->explim = CPFLOAT_EXPRANGE_STOR;
  }
}
free(xd);
free(intpartd);
free(refd);
free(oned);

float *xf = malloc(n * sizeof(*xf));
float *intpartf = malloc(n * sizeof(*intpartf));
float *reff = malloc(n * sizeof(*reff));
float *onef = malloc(n * sizeof(*onef));
for (size_t j = 0; j < n; j++)
  onef[j] = 1.0;

const int resf_isfinite [] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,  0, 0};
const int resf_isinf    [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1};
const int resf_isnan    [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0};
const int resf_isnormal [] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,  0, 0};
const int resf_fpclassify [] = {FP_NORMAL, FP_NORMAL,FP_NORMAL,
                                FP_NORMAL, FP_NORMAL, FP_NORMAL,
                                FP_NORMAL, FP_NORMAL, FP_NORMAL, FP_NORMAL,
                                FP_NAN, FP_SUBNORMAL, FP_SUBNORMAL,
                                FP_ZERO, FP_INFINITE, FP_INFINITE};

for (size_t i = 0; i < nformats; i++) {
  fpopts->precision = precision[i];
  fpopts->emax = emax[i];
  float xminf = minnormal(fpopts);
  float xminfs = minsubnormal(fpopts);
  float af[] = { -2., -1, -xminf, xminf, 1., CONST_SQRT2,
                CONST_LOG2E, 2., CONST_E, CONST_PI,
                nan_float(), -xminfs, xminfs, 0, -inf_float(), inf_float()};

  for (size_t i = 0; i < n_roundings; i++) {
    fpopts->round = det_roundings[i];

    cpfloatf(reff, af, n, fpopts);

    cpf_frexpf(xf, exp, af, n, fpopts);
    cpf_ldexpf(xf, xf, exp, n, fpopts);
    check_equality_float(xf, reff, n);

    cpf_modff(xf, intpartf, af, n, fpopts);
    for (size_t j = 0; j < n; j++)
      xf[j] = xf[j] + intpartf[j];
    check_equality_float(xf, reff, n);

    // These test would not work for +inf, -inf, and 0 because of the division.
    // Subnormals that would cause underflow are also removed.
    cpf_ilogbf(exp, af, n-5, fpopts);
    cpf_scalbnf(xf, onef, exp, n-5, fpopts);
    cpf_divf(xf, af, xf, n-5, fpopts);
    cpf_ldexpf(xf, xf, exp, n-5, fpopts);
    check_equality_float(xf, reff, n-5);

    cpf_ilogbf(exp, af, n-5, fpopts);
    for (size_t j = 0; j < n-5; j++)
      lexp[j] = (long)exp[j];
    cpf_scalblnf(xf, onef, lexp, n-5, fpopts);
    cpf_divf(xf, af, xf, n-5, fpopts);
    cpf_ldexpf(xf, xf, exp, n-5, fpopts);
    check_equality_float(xf, reff, n-5);

    // Removing nans and infs for this test.
    cpf_ilogbf(exp, af, n-6, fpopts);
    cpf_logbf(xf, af, n-6, fpopts);
    for (size_t j = 0; j < n-6; j++) {
      ck_assert(nan_safe_compare_float(xf[j], (float) exp[j]));
    }

    cpf_isfinitef(exp, af, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(both_true_or_false(exp[j], resf_isfinite[j]));

    cpf_isinff(exp, af, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(both_true_or_false(exp[j], resf_isinf[j]));

    cpf_isnanf(exp, af, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(both_true_or_false(exp[j], resf_isnan[j]));

    // Subnormals depend on the target format: we have to use the correct
    // exponent range and the correct subnormal handling mode to get the
    // reference results.
    fpopts->explim = CPFLOAT_EXPRANGE_TARG;
    fpopts->subnormal = CPFLOAT_SUBN_USE;
    cpf_isnormalf(exp, af, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(exp[j] == resf_isnormal[j]);
    cpf_fpclassifyf(exp, af, n, fpopts);
    for (size_t j = 0; j < n; j++)
      ck_assert(exp[j] == resf_fpclassify[j]);
    fpopts->subnormal = CPFLOAT_SUBN_RND;
    fpopts->explim = CPFLOAT_EXPRANGE_STOR;
  }
}
free(xf);
free(intpartf);
free(reff);
free(onef);
free(exp);
free(lexp);

#test enumeration_floating_point_numbers
printf("4c. Next and previous floating-point number\n");
fpopts->explim = CPFLOAT_EXPRANGE_TARG;
fpopts->subnormal = CPFLOAT_SUBN_USE;
// Subnormals
for (size_t mode=2; mode<3; mode++) {
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];

    size_t n = (ldexp(1.,fpopts->precision-1) - 1) + 2;
    double *ad = malloc(n * sizeof(*ad));
    double *xd = allocate_array_double(ad, n, mode);
    double *refd = malloc(n * sizeof(*refd));
    double *infd = malloc(n * sizeof(*infd));
    long double *infl = malloc(n * sizeof(*infl));
    for (size_t j = 0; j < n; j++) {
      infd[j] = inf_double();
      infl[j] = INFINITY;
    }

    init_fparray_double(refd, n, 0, minsubnormal(fpopts)); // 0 -> min. normal.
    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad, infd, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad, infl, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);
    csign_intarray_double((uint64_t *)infd, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad+1, infd, n-1, fpopts);
    check_equality_double(xd, refd, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad+1, infl, n-1, fpopts);
    check_equality_double(xd, refd, n-1);

    init_fparray_double(refd, n, 0, minsubnormal(fpopts)); // 0 -> min. normal.
    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad+1, infd, n-1, fpopts);
    check_equality_double(xd, refd, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad+1, infl, n-1, fpopts);
    check_equality_double(xd, refd, n-1);
    csign_intarray_double((uint64_t *) infd, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad, infd, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad, infl, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);

    free(ad);
    free(refd);
    free_array_double(xd, mode);
    free(infd);

    float *af = malloc(n * sizeof(* af));
    float *xf = allocate_array_float(af, n, mode);
    float *reff = malloc(n * sizeof(* reff));
    float *inff = malloc(n * sizeof(*inff));
    for (size_t j = 0; j < n; j++)
      inff[j] = inf_float();

    init_fparray_float(reff, n, 0, minsubnormal(fpopts)); // 0 -> min. normal.
    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af, inff, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af, infl, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);
    csign_intarray_float((uint32_t *)inff, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af+1, inff, n-1, fpopts);
    check_equality_float(xf, reff, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af+1, infl, n-1, fpopts);
    check_equality_float(xf, reff, n-1);

    init_fparray_float(reff, n, 0, minsubnormal(fpopts)); // -min. normal. -> 0
    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af+1, inff, n-1, fpopts);
    check_equality_float(xf, reff, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af+1, infl, n-1, fpopts);
    check_equality_float(xf, reff, n-1);
    csign_intarray_float((uint32_t *) inff, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af, inff, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af, infl, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);

    free(af);
    free_array_float(xf, mode);
    free(inff);
    free(infl);
  }
}

// Normals
for (size_t mode=2; mode<3; mode++) {
  for (size_t i = 0; i < nformats; i++) {
    fpopts->precision = precision[i];
    fpopts->emax = emax[i];

    size_t n = ldexp(1., fpopts->precision-1) * 2 * fpopts->emax + 2;
    double *ad = malloc(n * sizeof(*ad));
    double *xd = allocate_array_double(ad, n, mode);
    double *refd = malloc(n * sizeof(*refd));
    double *infd = malloc(n * sizeof(*infd));
    long double *infl = malloc(n * sizeof(*infl));
    for (size_t j = 0; j < n; j++) {
      infd[j] = inf_double();
      infl[j] = INFINITY;
    }

    init_intarray_double((uint64_t *)(refd+1), n-1, intminnormal_double(fpopts),
                         1ul << (52 - fpopts->precision + 1));
    refd[0] = maxsubnormal(fpopts);
    refd[n-1] = inf_double();

    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad, infd, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad, infl, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);
    csign_intarray_double((uint64_t *)infd, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad+1, infd, n-1, fpopts);
    check_equality_double(xd, refd, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad+1, infl, n-1, fpopts);
    check_equality_double(xd, refd, n-1);

    csign_intarray_double((uint64_t *)refd, n);
    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad, infd, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad, infl, n-1, fpopts);
    check_equality_double(xd, refd+1, n-1);
    csign_intarray_double((uint64_t *)infd, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_double(ad, refd, n);
    cpf_nextafter(xd, ad+1, infd, n-1, fpopts);
    check_equality_double(xd, refd, n-1);
    copy_array_double(ad, refd, n);
    cpf_nexttoward(xd, ad+1, infl, n-1, fpopts);
    check_equality_double(xd, refd, n-1);

    free(ad);
    free(refd);
    free_array_double(xd, mode);
    free(infd);

    float *af = malloc(n * sizeof(* af));
    float *xf = allocate_array_float(af, n, mode);
    float *reff = malloc(n * sizeof(* reff));
    float *inff = malloc(n * sizeof(*inff));
    for (size_t j = 0; j < n; j++)
      inff[j] = inf_float();

    init_intarray_float((uint32_t *)(reff+1), n-1, intminnormal_float(fpopts),
                         1ul << (23 - fpopts->precision + 1));
    reff[0] = maxsubnormal(fpopts);
    reff[n-1] = inf_double();

    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af, inff, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af, infl, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);
    csign_intarray_float((uint32_t *)inff, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af+1, inff, n-1, fpopts);
    check_equality_float(xf, reff, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af+1, infl, n-1, fpopts);
    check_equality_float(xf, reff, n-1);

    csign_intarray_float((uint32_t *) reff, n);
    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af, inff, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af, infl, n-1, fpopts);
    check_equality_float(xf, reff+1, n-1);
    csign_intarray_float((uint32_t *)inff, n);
    for (size_t j = 0; j < n; j++) infl[j] = -infl[j];
    copy_array_float(af, reff, n);
    cpf_nextafterf(xf, af+1, inff, n-1, fpopts);
    check_equality_float(xf, reff, n-1);
    copy_array_float(af, reff, n);
    cpf_nexttowardf(xf, af+1, infl, n-1, fpopts);
    check_equality_float(xf, reff, n-1);

    free(af);
    free_array_float(xf, mode);
    free(inff);
    free(infl);
  }
}

#test integer_rounding
printf("4c. Integer rounding\n");

fpopts->emax = emax[1];
fpopts->precision = precision[1];

size_t n = 32;
double *xd = malloc(n * sizeof(*xd));
int *xi = malloc(n * sizeof(*xi));
long *xl = malloc(n * sizeof(*xl));
long long *xll = malloc(n * sizeof(*xll));
double ad[]  = {-3.00, -2.75, -2.50, -2.25,
                -2.00, -1.75, -1.50, -1.25,
                -1.00, -0.75, -0.50, -0.25,
                -0.00,
                0.00,  0.25,  0.50,  0.75,
                1.00,  1.25,  1.50,  1.75,
                2.00,  2.25,  2.50,  2.75,
                -1234.56789123456789,
                -1234.36789123456789,
                1234.56789123456789,
                1234.36789123456789,
                -INFINITY, INFINITY, NAN};

fpopts->round = CPFLOAT_RND_NA;
double rd1[] = {-3.00, -3.00, -3.00, -2.00,
                -2.00, -2.00, -2.00, -1.00,
                -1.00, -1.00, -1.00, -0.00,
                -0.00,
                0.00,  0.00,  1.00,  1.00,
                1.00,  1.00,  2.00,  2.00,
                2.00,  2.00,  3.00,  3.00,
                -1235.00, -1234.00, 1235.00, 1234.00,
                -INFINITY, INFINITY, NAN};

cpf_round(xd, ad, n, fpopts);
check_equality_double(rd1, xd, n);
cpf_lround(xl, ad, n-3, fpopts);
check_equality_double_long(rd1, xl, n-3);
cpf_llround(xll, ad, n-3, fpopts);
check_equality_double_long_long(rd1, xll, n-3);

cpf_rint(xd, xi, ad, n, fpopts);
check_equality_double(rd1, xd, n);
cpf_nearbyint(xd, ad, n, fpopts);
check_equality_double(rd1, xd, n);
cpf_lrint(xl, xi, ad, n-3, fpopts);
check_equality_double_long(rd1, xl, n-3);
cpf_llrint(xll, xi, ad, n-3, fpopts);
check_equality_double_long_long(rd1, xll, n-3);

fpopts->round = CPFLOAT_RND_NZ;
double rd2[] = {-3.00, -3.00, -2.00, -2.00,
                -2.00, -2.00, -1.00, -1.00,
                -1.00, -1.00, -0.00, -0.00,
                -0.00,
                0.00,  0.00,  0.00,  1.00,
                1.00,  1.00,  1.00,  2.00,
                2.00,  2.00,  2.00,  3.00,
                -1235.00, -1234.00, 1235.00, 1234.00,
                -INFINITY, INFINITY, NAN};
cpf_rint(xd, xi, ad, n, fpopts);
check_equality_double(rd2, xd, n);
cpf_nearbyint(xd, ad, n, fpopts);
check_equality_double(rd2, xd, n);
cpf_lrint(xl, xi, ad, n-3, fpopts);
check_equality_double_long(rd2, xl, n-3);
cpf_llrint(xll, xi, ad, n-3, fpopts);
check_equality_double_long_long(rd2, xll, n-3);

fpopts->round = CPFLOAT_RND_NE;
double rd3[] = {-3.00, -3.00, -2.00, -2.00,
                -2.00, -2.00, -2.00, -1.00,
                -1.00, -1.00, -0.00, -0.00,
                -0.00,
                0.00,  0.00,  0.00,  1.00,
                1.00,  1.00,  2.00,  2.00,
                2.00,  2.00,  2.00,  3.00,
                -1235.00, -1234.00, 1235.00, 1234.00,
                -INFINITY, INFINITY, NAN};
cpf_rint(xd, xi, ad, n, fpopts);
check_equality_double(rd3, xd, n);
cpf_nearbyint(xd, ad, n, fpopts);
check_equality_double(rd3, xd, n);
cpf_lrint(xl, xi, ad, n-3, fpopts);
check_equality_double_long(rd3, xl, n-3);
cpf_llrint(xll, xi, ad, n-3, fpopts);
check_equality_double_long_long(rd3, xll, n-3);

fpopts->round = CPFLOAT_RND_TP;
double rd4[] = {-3.00, -2.00, -2.00, -2.00,
                -2.00, -1.00, -1.00, -1.00,
                -1.00, -0.00, -0.00, -0.00,
                -0.00,
                0.00,  1.00,  1.00,  1.00,
                1.00,  2.00,  2.00,  2.00,
                2.00,  3.00,  3.00,  3.00,
                -1234.00, -1234.00, 1235.00, 1235.00,
                -INFINITY, INFINITY, NAN};
cpf_rint(xd, xi, ad, n, fpopts);
check_equality_double(rd4, xd, n);
cpf_nearbyint(xd, ad, n, fpopts);
check_equality_double(rd4, xd, n);
cpf_ceil(xd, ad, n, fpopts);
check_equality_double(rd4, xd, n);
cpf_lrint(xl, xi, ad, n-3, fpopts);
check_equality_double_long(rd4, xl, n-3);
cpf_llrint(xll, xi, ad, n-3, fpopts);
check_equality_double_long_long(rd4, xll, n-3);

fpopts->round = CPFLOAT_RND_TN;
double rd5[] = {-3.00, -3.00, -3.00, -3.00,
                -2.00, -2.00, -2.00, -2.00,
                -1.00, -1.00, -1.00, -1.00,
                -0.00,
                0.00,  0.00,  0.00,  0.00,
                1.00,  1.00,  1.00,  1.00,
                2.00,  2.00,  2.00,  2.00,
                -1235.00, -1235.00, 1234.00, 1234.00,
                -INFINITY, INFINITY, NAN};
cpf_rint(xd, xi, ad, n, fpopts);
check_equality_double(rd5, xd, n);
cpf_nearbyint(xd, ad, n, fpopts);
check_equality_double(rd5, xd, n);
cpf_floor(xd, ad, n, fpopts);
check_equality_double(rd5, xd, n);
cpf_lrint(xl, xi, ad, n-3, fpopts);
check_equality_double_long(rd5, xl, n-3);
cpf_llrint(xll, xi, ad, n-3, fpopts);
check_equality_double_long_long(rd5, xll, n-3);

fpopts->round = CPFLOAT_RND_TZ;
double rd6[] = {-3.00, -2.00, -2.00, -2.00,
                -2.00, -1.00, -1.00, -1.00,
                -1.00, -0.00, -0.00, -0.00,
                -0.00,
                0.00,  0.00,  0.00,  0.00,
                1.00,  1.00,  1.00,  1.00,
                2.00,  2.00,  2.00,  2.00,
                -1234.00, -1234.00, 1234.00, 1234.00,
                -INFINITY, INFINITY, NAN};
cpf_rint(xd, xi, ad, n, fpopts);
check_equality_double(rd6, xd, n);
cpf_nearbyint(xd, ad, n, fpopts);
check_equality_double(rd6, xd, n);
cpf_trunc(xd, ad, n, fpopts);
check_equality_double(rd6, xd, n);
cpf_lrint(xl, xi, ad, n-3, fpopts);
check_equality_double_long(rd6, xl, n-3);
cpf_llrint(xll, xi, ad, n-3, fpopts);
check_equality_double_long_long(rd6, xll, n-3);

fpopts->round = CPFLOAT_RND_OD;
double rd7[] = {-3.00, -3.00, -3.00, -3.00,
                -2.00, -1.00, -1.00, -1.00,
                -1.00, -1.00, -1.00, -1.00,
                -0.00,
                0.00,  1.00,  1.00,  1.00,
                1.00,  1.00,  1.00,  1.00,
                2.00,  3.00,  3.00,  3.00,
                -1235.00, -1235.00, 1235.00, 1235.00,
                -INFINITY, INFINITY, NAN};
cpf_rint(xd, xi, ad, n, fpopts);
check_equality_double(rd7, xd, n);
cpf_nearbyint(xd, ad, n, fpopts);
check_equality_double(rd7, xd, n);
cpf_lrint(xl, xi, ad, n-3, fpopts);
check_equality_double_long(rd7, xl, n-3);
cpf_llrint(xll, xi, ad, n-3, fpopts);
check_equality_double_long_long(rd7, xll, n-3);


#test arithmetic_operations_large
printf("4d. Arithmetic operations on large arrays\n");
size_t i = 0;
fpopts->precision = precision[i];
fpopts->emax = emax[i];
size_t n = 3 * (ldexp(1., fpopts->precision-1) * (ldexp(1., 5) - 2) - 1) - 2;
double *ad = malloc(n * sizeof(*ad));
double *bd = malloc(n * sizeof(*bd));
double *cd = malloc(n * sizeof(*cd));
double *refd = malloc(n * sizeof(*refd));
double *xd = malloc(n * sizeof(*xd));
uint64_t stepd = 1ul << (52-fpopts->precision + 1);
init_intarray_rounding_double((uint64_t *)ad, n,
                              intminnormal_double(fpopts), stepd);
init_intarray_rounding_double((uint64_t *)bd, n,
                              intminnormal_double(fpopts)+stepd, stepd);
init_intarray_rounding_double((uint64_t *)cd, n,
                              intminnormal_double(fpopts)+2*stepd, stepd);
test_univariate_arith_op_double(xd, refd, ad, n,
                                fpopts, cpf_cos, cos);
test_bivariate_arith_op_double(xd, refd, ad, bd, n,
                               fpopts, cpf_add, ref_add);
test_trivariate_arith_op_double(xd, refd, ad, bd, cd, n,
                                fpopts, cpf_fma, fma);
free(ad);
free(bd);
free(cd);
free(refd);
free(xd);

float *af = malloc(n * sizeof(*af));
float *bf = malloc(n * sizeof(*bf));
float *cf = malloc(n * sizeof(*cf));
float *reff = malloc(n * sizeof(*reff));
float *xf = malloc(n * sizeof(*xf));
uint32_t stepf = 1ul << (52-fpopts->precision + 1);
init_intarray_rounding_float((uint32_t *)af, n,
                             intminnormal_float(fpopts), stepf);
init_intarray_rounding_float((uint32_t *)bf, n,
                             intminnormal_float(fpopts)+stepf, stepf);
init_intarray_rounding_float((uint32_t *)cf, n,
                             intminnormal_float(fpopts)+2*stepf, stepf);
test_univariate_arith_op_float(xf, reff, af, n,
                               fpopts, cpf_cosf, cosf);
test_bivariate_arith_op_float(xf, reff, af, bf, n,
                              fpopts, cpf_addf, ref_addf);
test_trivariate_arith_op_float(xf, reff, af, bf, cf, n,
                               fpopts, cpf_fmaf, fmaf);
free(af);
free(bf);
free(cf);
free(reff);
free(xf);

#main-pre
putenv("CK_VERBOSITY=verbose");
tcase_add_unchecked_fixture(tc1_1, fpopts_setup, fpopts_teardown);
tcase_add_unchecked_fixture(tc1_2, fpopts_setup, fpopts_teardown);
tcase_add_unchecked_fixture(tc1_3, fpopts_setup, fpopts_teardown);
tcase_add_unchecked_fixture(tc1_4, fpopts_setup, fpopts_teardown);
tcase_set_timeout(tc1_1, 1800);
tcase_set_timeout(tc1_2, 1800);
tcase_set_timeout(tc1_3, 1800);
tcase_set_timeout(tc1_4, 1800);
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
