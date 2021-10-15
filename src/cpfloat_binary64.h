/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/**
 * @file cpfloat_binary64.h
 * @brief CPFloat functions for `double` arrays.
 */

#ifndef _CPFLOAT_BINARY64_
#define _CPFLOAT_BINARY64_

#include "cpfloat_definitions.h"
#include "cpfloat_docmacros.h"

/* Validation of floating-point parameters. */
doc_cpfloat_validate_optstruct(double, 26, 53, 1023)
static inline int cpfloat_validate_optstruct(const optstruct *fpopts);

/* Rounding functions. */
doc_cpfloat(double, 53, 1023)
static inline int cpfloat(double *X, const double *A, const size_t numelem,
                          optstruct *fpopts);
doc_cpfloat(double, 53, 1023)
static inline int cpf_fpround(double *X, const double *A,
                                  const size_t numelem, optstruct *fpopts);

/* Elementary arithmetic operations. */
doc_cpf_bivariate(sum, X[i] = A[i] + B[i], 53, 1023)
static inline int cpf_add(double *X, const double *A, const double *B,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(difference, X[i] = A[i] - B[i], 53, 1023)
static inline int cpf_sub(double *X, const double *A, const double *B,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(product, X[i] = A[i] * B[i], 53, 1023)
static inline int cpf_mul(double *X, const double *A, const double *B,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(ratio, X[i] = A[i] / B[i], 53, 1023)
static inline int cpf_div(double *X, const double *A, const double *B,
                              const size_t numelem, optstruct *fpopts);

/* Trigonometric functions. */
doc_cpf_univariate(trigonometric cosine, X[i] = cos(A[i]), 53, 1023)
static inline int cpf_cos(double *X, const double *A,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(trigonometric sine, X[i] = sin(A[i]), 53, 1023)
static inline int cpf_sin(double *X, const double *A,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(trigonometric tangent, X[i] = tan(A[i]), 53, 1023)
static inline int cpf_tan(double *X, const double *A,
                              const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(inverse trigonometric cosine,
                       X[i] = acos(A[i]), 53, 1023)
static inline int cpf_acos(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse trigonometric sine,
                       X[i] = asin(A[i]), 53, 1023)
static inline int cpf_asin(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse trigonometric tangent,
                       X[i] = atan(A[i]), 53, 1023)
static inline int cpf_atan(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(2-argument arctangent,
                      X[i] = atan2 (A[i], B[i]), 53, 1023)
static inline int cpf_atan2(double *X, const double *A, const double *B,
                                const size_t numelem, optstruct *fpopts);

/* Hyperbolic functions. */
doc_cpf_univariate(hyperbolic cosine, X[i] = cosh(A[i]), 53, 1023)
static inline int cpf_cosh(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(hyperbolic sine, X[i] = sinh(A[i]), 53, 1023)
static inline int cpf_sinh(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(hyperbolic tangent , X[i] = tanh(A[i]), 53, 1023)
static inline int cpf_tanh(double *X, const double *A,
                              const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(inverse hyperbolic cosine,
                       X[i] = arcosh(A[i]), 53, 1023)
static inline int cpf_acosh(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse hyperbolic sine,
                       X[i] = arsinh(A[i]), 53, 1023)
static inline int cpf_asinh(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse hyperbolic tangent,
                       X[i] = artanh(A[i]), 53, 1023)
static inline int cpf_atanh(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);

/* Exponentiation and logarithmic functions. */
doc_cpf_univariate(exponential, X[i] = exp(A[i]), 53, 1023)
static inline int cpf_exp(double *X, const double *A,
                              const size_t numelem, optstruct *fpopts);

doc_cpf_frexp(53, 1023)
static inline int cpf_frexp(double *X, int *exp, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(2, 53, 1023)
static inline int cpf_ldexp(double *X, const double *A, const int *exp,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(natural logarithm, X[i] = log(A[i]), 53, 1023)
static inline int cpf_log(double *X, const double *A,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-10 logarithm, X[i] = log10(A[i]), 53, 1023)
static inline int cpf_log10(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_modf(53, 1023)
static inline int cpf_modf(double *X, double *intpart, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-2 exponential, X[i] = 2^(A[i]), 53, 1023)
static inline int cpf_exp2(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(exp(x) - 1, X[i] = exp(A[i]) - 1, 53, 1023)
static inline int cpf_expm1(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_ilogb(53, 1023)
static inline int cpf_ilogb(int *exp, const double *A,
                                const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(natural logarithm of number shifted by one,
                       X[i] = log(1+A[i]), 53, 1023)
static inline int cpf_log1p(double *X, const double *A,
                                size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-2 logarithm, X[i] = log2(A[i]), 53, 1023)
static inline int cpf_log2(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(FLT_RADIX, 53, 1023)
  static inline int cpf_scalbn(double *X, const double *A, const int *exp,
                                   const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(FLT_RADIX, 53, 1023)
  static inline int cpf_scalbln(double *X, const double *A,
                                    const long int *exp, const size_t numelem,
                                    optstruct *fpopts);

/* Power functions. */
doc_cpf_bivariate(real powers, X[i] = A[i]^B[i], 53, 1023)
static inline int cpf_pow(double *X, const double *A, const double *B,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(square root, X[i] = sqrt(A[i]), 53, 1023)
static inline int cpf_sqrt(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(cube root, X[i] = cbrt(A[i]), 53, 1023)
static inline int cpf_cbrt(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(hypotenuse of a right-angle triangle,
                      X[i] = sqrt(A[i]^2 + B[i]^2), 53, 1023)
static inline int cpf_hypot(double *X, const double *A, const double *B,
                                const size_t numelem, optstruct *fpopts);

/* Error and gamma functions. */
doc_cpf_univariate(error function, X[i] = erf(A[i]), 53, 1023)
static inline int cpf_erf(double *X, const double *A,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(complementary error function, X[i] = erfc(A[i]), 53, 1023)
static inline int cpf_erfc(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(gamma function, X[i] = gamma(A[i]), 53, 1023)
static inline int cpf_tgamma(double *X, const double *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(natural logarithm of absolute value of gamma function,
                       X[i] = log(abs(gamma(A[i]))), 53, 1023)
static inline int cpf_lgamma(double *X, const double *A,
                                 const size_t numelem, optstruct *fpopts);

/* Rounding and remainder functions. */
doc_cpf_univariate(ceiling function, X[i] = ceil(A[i]), 53, 1023)
static inline int cpf_ceil(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(floor function, X[i] = floor(A[i]), 53, 1023)
static inline int cpf_floor(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(integer truncation, X[i] = trunc(A[i]), 53, 1023)
static inline int cpf_trunc(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(closest integer (with round-to-nearest),
                       X[i] = round(A[i]), 53, 1023)
static inline int cpf_round(double *X, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(closest integer (with round-to-nearest),
                       X[i] = round(A[i]), 53, 1023)
static inline int cpf_lround(long *X, const double *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_univariate_nobitflip(closest integer (with round-to-nearest),
                       X[i] = round(A[i]), 53, 1023)
static inline int cpf_llround(long long *X, const double *A,
                                 const size_t numelem, optstruct *fpopts);

doc_cpf_rint(PMAX, EMAX)
static inline int cpf_rint(double *X, int *exception, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_rint(PMAX, EMAX)
static inline int cpf_lrint(long *X, int *exception, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_rint(PMAX, EMAX)
static inline int cpf_llrint(long long *X, int *exception, const double *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_nearbyint(PMAX, EMAX)
static inline int cpf_nearbyint(double *X, const double *A,
                                    const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(remainder of the floating point division,
                      X[i] = A[i]^2 - k * B[i]
                      for largest k such that k * B[i] < A[i], 53, 1023)
static inline int cpf_remainder(double *X, const double *A, const double *B,
                                    const size_t numelem, optstruct *fpopts);

doc_cpf_remquo(PMAX, EMAX)
static inline int cpf_remquo(double *X, int *quot,
                                 const double *A, const double *B,
                                 const size_t numelem, optstruct *fpopts);

/* Floating-point manipulation functions. */
doc_cpf_bivariate(number from magnitude and sign,
                      X[i] = sign(A[i]) * abs(B[i]), 53, 1023)
static inline int cpf_copysign(double *X, const double *A, const double *B,
                                   const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(next floating-point number in specified direction,
                      the floating-point number closest to A[i] in the
                      direction of B[i], 53, 1023)
static inline int cpf_nextafter(double *X, double *A, const double *B,
                                    const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(next floating-point number in specified direction,
                      the floating-point number closest to A[i] in the
                      direction of B[i], 53, 1023)
static inline int cpf_nexttoward(double *X, double *A,
                                     const long double *B, const size_t numelem,
                                     optstruct *fpopts);

/* Minimum, maximum, difference functions. */
doc_cpf_bivariate(positive difference, X[i] = abs(A[i]) - B[i], 53, 1023)
static inline int cpf_fdim(double *X, const double *A, const double *B,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(element-wise maximum, X[i] = max(A[i], B[i]), 53, 1023)
static inline int cpf_fmax(double *X, const double *A, const double *B,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(element-wise minimum, X[i] = min(A[i], B[i]), 53, 1023)
static inline int cpf_fmin(double *X, const double *A, const double *B,
                               const size_t numelem, optstruct *fpopts);

/* Classification. */
doc_cpf_fpclassify(53, 1023)
static inline int cpf_fpclassify(int *r, const double *A,
                                     const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(finite, 53, 1023)
static inline int cpf_isfinite(int *r, const double *A,
                                   const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(infinite, 53, 1023)
static inline int cpf_isinf(int *r, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(not a number, 53, 1023)
static inline int cpf_isnan(int *r, const double *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(normal, 53, 1023)
static inline int cpf_isnormal(int *r, const double *A,
                                   const size_t numelem, optstruct *fpopts);

/* Other functions. */
doc_cpf_univariate(absolute value, X[i] = abs(A[i]), 53, 1023)
static inline int cpf_fabs(double *X, const double *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_trivariate(fused multiply-add , X[i] = A[i] * B[i] + C[i], 53, 1023)
static inline int cpf_fma(double *X, const double *A, const double *B,
                              const double *C, const size_t numelem,
                              optstruct *fpopts);

/** @cond */
#define FUNSUFFIX
#define FPTYPE double
#define INTTYPE uint64_t
#define INTSUFFIX ULL
#define DEFPREC 53
#define DEFEMAX 1023
#define DEFEMIN -1022
#define NLEADBITS 12
#define NBITS 64
#define FULLMASK 0xFFFFFFFFFFFFFFFFULL
#define ABSMASK 0x7FFFFFFFFFFFFFFFULL
#define SIGNMASK 0x8000000000000000ULL
#define EXPMASK 0x7FF0000000000000ULL
#define FRACMASK 0x000FFFFFFFFFFFFFULL

#ifdef PCG_VARIANTS_H_INCLUDED
#define INITRAND(seed) pcg64_srandom_r(seed, time(NULL), (intptr_t)seed);
#define ADVANCERAND(seed, thread, nloc) pcg64_advance_r(seed, thread *nloc - 1);
#define GENRAND(seed) pcg64_random_r(seed)
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#warning "The default C random number generator is being used."
#warning "Please compile with -include <path-to-pcg_variants.h>"
#warning "and link with -L <path-to-libpcg_random.a> -lpcg_random."
#define INITRAND(seed) *seed = time(NULL);
#define GEN_SINGLE_RAND(seed)                                                  \
  ((INTTYPE)rand_r((unsigned int *)seed) +                                     \
   ((INTTYPE)rand_r((unsigned int *)seed) << 31))
#endif /* #ifdef PCG_VARIANTS_H_INCLUDED */

#include "cpfloat_threshold_binary64.h"
#include "cpfloat_template.h"
/** @endcond */

#endif  /* #ifndef _CPFLOAT_BINARY64_ */

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
