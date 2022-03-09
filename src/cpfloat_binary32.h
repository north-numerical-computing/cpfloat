/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/**
 * @file cpfloat_binary32.h
 * @brief CPFloat functions for `float` arrays.
 */

#ifndef _CPFLOAT_BINARY32_
#define _CPFLOAT_BINARY32_

#include "cpfloat_definitions.h"
#include "cpfloat_docmacros.h"

/* Validation of floating-point parameters. */
doc_cpfloat_validate_optstruct(double, 12, 24, 127)
static inline int cpfloat_validate_optstructf(const optstruct *fpopts);

/* Rounding functions. */
doc_cpfloat(float, 24, 127)
static inline int cpfloatf(float *X, const float *A, const size_t numelem,
                           optstruct *fpopts);
doc_cpfloat(float, 24, 127)
static inline int cpf_fproundf(float *X, const float *A,
                                   const size_t numelem, optstruct *fpopts);

/* Elementary arithmetic operations. */
doc_cpf_bivariate(sum, X[i] = A[i] + B[i], 24, 127)
static inline int cpf_addf(float *X, const float *A, const float *B,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(difference, X[i] = A[i] - B[i], 24, 127)
static inline int cpf_subf(float *X, const float *A, const float *B,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(product, X[i] = A[i] * B[i], 24, 127)
static inline int cpf_mulf(float *X, const float *A, const float *B,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(ratio, X[i] = A[i] / B[i], 24, 127)
static inline int cpf_divf(float *X, const float *A, const float *B,
                               const size_t numelem, optstruct *fpopts);

/* Trigonometric functions. */
doc_cpf_univariate(trigonometric cosine, X[i] = cos(A[i]), 24, 127)
static inline int cpf_cosf(float *X, const float *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(trigonometric sine, X[i] = sin(A[i]), 24, 127)
static inline int cpf_sinf(float *X, const float *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(trigonometric tangent, X[i] = tan(A[i]), 24, 127)
static inline int cpf_tanf(float *X, const float *A,
                               const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(inverse trigonometric cosine,
                       X[i] = acos(A[i]), 24, 127)
static inline int cpf_acosf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse trigonometric sine,
                       X[i] = asin(A[i]), 24, 127)
static inline int cpf_asinf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse trigonometric tangent,
                       X[i] = atan(A[i]), 24, 127)
static inline int cpf_atanf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(2-argument arctangent,
                      X[i] = atan2 (A[i], B[i]), 24, 127)
static inline int cpf_atan2f(float *X, const float *A, const float *B,
                                 const size_t numelem, optstruct *fpopts);

/* Hyperbolic functions. */
doc_cpf_univariate(hyperbolic cosine, X[i] = cosh(A[i]), 24, 127)
static inline int cpf_coshf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(hyperbolic sine, X[i] = sinh(A[i]), 24, 127)
static inline int cpf_sinhf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(hyperbolic tangent , X[i] = tanh(A[i]), 24, 127)
static inline int cpf_tanhf(float *X, const float *A,
                               const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(inverse hyperbolic cosine,
                       X[i] = arcosh(A[i]), 24, 127)
static inline int cpf_acoshf(float *X, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse hyperbolic sine,
                       X[i] = arsinh(A[i]), 24, 127)
static inline int cpf_asinhf(float *X, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse hyperbolic tangent,
                       X[i] = artanh(A[i]), 24, 127)
static inline int cpf_atanhf(float *X, const float *A,
                                 const size_t numelem, optstruct *fpopts);

/* Exponentiation and logarithmic functions. */
doc_cpf_univariate(exponential, X[i] = exp(A[i]), 24, 127)
static inline int cpf_expf(float *X, const float *A,
                               const size_t numelem, optstruct *fpopts);

doc_cpf_frexp(24, 127)
static inline int cpf_frexpf(float *X, int *exp, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(2, 24, 127)
static inline int cpf_ldexpf(float *X, const float *A, const int *exp,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(natural logarithm, X[i] = log(A[i]), 24, 127)
static inline int cpf_logf(float *X, const float *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-10 logarithm, X[i] = log10(A[i]), 24, 127)
static inline int cpf_log10f(float *X, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_modf(24, 127)
static inline int cpf_modff(float *X, float *intpart, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-2 exponential, X[i] = 2^(A[i]), 24, 127)
static inline int cpf_exp2f(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(exp(x) - 1, X[i] = exp(A[i]) - 1, 24, 127)
static inline int cpf_expm1f(float *X, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_ilogb(24, 127)
static inline int cpf_ilogbf(int *exp, const float *A,
                                 const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(natural logarithm of number shifted by one,
                       X[i] = log(1+A[i]), 24, 127)
static inline int cpf_log1pf(float *X, const float *A,
                                 size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-2 logarithm, X[i] = log2(A[i]), 24, 127)
static inline int cpf_log2f(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-FLT_RADIX logarithm of absolute value,
                   X[i] = log(abs(A[i])), 24, 127)
static inline int cpf_logbf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(FLT_RADIX, 24, 127)
  static inline int cpf_scalbnf(float *X, const float *A, const int *exp,
                                    const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(FLT_RADIX, 24, 127)
  static inline int cpf_scalblnf(float *X, const float *A,
                                     const long int *exp, const size_t numelem,
                                     optstruct *fpopts);

/* Power functions. */
doc_cpf_bivariate(real powers, X[i] = A[i]^B[i], 24, 127)
static inline int cpf_powf(float *X, const float *A, const float *B,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(square root, X[i] = sqrt(A[i]), 24, 127)
static inline int cpf_sqrtf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(cube root, X[i] = cbrt(A[i]), 24, 127)
static inline int cpf_cbrtf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(hypotenuse of a right-angle triangle,
                      X[i] = sqrt(A[i]^2 + B[i]^2), 24, 127)
static inline int cpf_hypotf(float *X, const float *A, const float *B,
                                 const size_t numelem, optstruct *fpopts);

/* Error and gamma functions. */
doc_cpf_univariate(error function, X[i] = erf(A[i]), 24, 127)
static inline int cpf_erff(float *X, const float *A,
                               const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(complementary error function, X[i] = erfc(A[i]), 24, 127)
static inline int cpf_erfcf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(gamma function, X[i] = gamma(A[i]), 24, 127)
static inline int cpf_tgammaf(float *X, const float *A,
                                  const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(natural logarithm of absolute value of gamma function,
                       X[i] = log(abs(gamma(A[i]))), 24, 127)
static inline int cpf_lgammaf(float *X, const float *A,
                                  const size_t numelem, optstruct *fpopts);

/* Rounding and remainder functions. */
doc_cpf_univariate(ceiling function, X[i] = ceil(A[i]), 24, 127)
static inline int cpf_ceilf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(floor function, X[i] = floor(A[i]), 24, 127)
static inline int cpf_floorf(float *X, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(floating-point remainder of division,
                  X[i] = A[i] mod B[i], 24, 127)
static inline int cpf_fmodf(float *X, const float *A, const float *B,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(integer truncation, X[i] = trunc(A[i]), 24, 127)
static inline int cpf_truncf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(closest integer (with round-to-nearest),
                       X[i] = round(A[i]), 24, 127)
static inline int cpf_roundf(float *X, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(closest integer (with round-to-nearest),
                       X[i] = round(A[i]), 24, 127)
static inline int cpf_lroundf(long *X, const float *A,
                                  const size_t numelem, optstruct *fpopts);
doc_cpf_univariate_nobitflip(closest integer (with round-to-nearest),
                       X[i] = round(A[i]), 24, 127)
static inline int cpf_llroundf(long long *X, const float *A,
                                   const size_t numelem, optstruct *fpopts);

doc_cpf_rint(PMAX, EMAX)
static inline int cpf_rintf(float *X, int *exception, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_rint(PMAX, EMAX)
static inline int cpf_lrintf(long *X, int *exception, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_rint(PMAX, EMAX)
static inline int cpf_llrintf(long long *X, int *exception, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_nearbyint(PMAX, EMAX)
static inline int cpf_nearbyintf(float *X, const float *A,
                                     const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(remainder of the floating point division,
                      X[i] = A[i]^2 - k * B[i]
                      for largest k such that k * B[i] < A[i], 24, 127)
static inline int cpf_remainderf(float *X, const float *A, const float *B,
                                     const size_t numelem, optstruct *fpopts);

doc_cpf_remquo(PMAX, EMAX)
static inline int cpf_remquof(float *X, int *quot,
                                  const float *A, const float *B,
                                  const size_t numelem, optstruct *fpopts);

/* Floating-point manipulation functions. */
doc_cpf_bivariate(number from magnitude and sign,
                      X[i] = sign(A[i]) * abs(B[i]), 24, 127)
static inline int cpf_copysignf(float *X, const float *A, const float *B,
                                    const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(next floating-point number in specified direction,
                      the floating-point number closest to A[i] in the
                      direction of B[i], 24, 127)
static inline int cpf_nextafterf(float *X, const float *A, const float *B,
                                     const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(next floating-point number in specified direction,
                      the floating-point number closest to A[i] in the
                      direction of B[i], 24, 127)
static inline int cpf_nexttowardf(float *X, const float *A,
                                      const long double *B,
                                      const size_t numelem,
                                      optstruct *fpopts);

/* Minimum, maximum, difference functions. */
doc_cpf_bivariate(positive difference, X[i] = abs(A[i]) - B[i], 24, 127)
static inline int cpf_fdimf(float *X, const float *A, const float *B,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(element-wise maximum, X[i] = max(A[i], B[i]), 24, 127)
static inline int cpf_fmaxf(float *X, const float *A, const float *B,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(element-wise minimum, X[i] = min(A[i], B[i]), 24, 127)
static inline int cpf_fminf(float *X, const float *A, const float *B,
                                const size_t numelem, optstruct *fpopts);

/* Classification. */
doc_cpf_fpclassify(24, 127)
static inline int cpf_fpclassifyf(int *r, const float *A,
                                      const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(finite, 24, 127)
static inline int cpf_isfinitef(int *r, const float *A,
                                    const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(infinite, 24, 127)
static inline int cpf_isinff(int *r, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(not a number, 24, 127)
static inline int cpf_isnanf(int *r, const float *A,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_isfun(normal, 24, 127)
static inline int cpf_isnormalf(int *r, const float *A,
                                    const size_t numelem, optstruct *fpopts);

/* Other functions. */
doc_cpf_univariate(absolute value, X[i] = abs(A[i]), 24, 127)
static inline int cpf_fabsf(float *X, const float *A,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_trivariate(fused multiply-add , X[i] = A[i] * B[i] + C[i], 24, 127)
static inline int cpf_fmaf(float *X, const float *A, const float *B,
                               const float *C, const size_t numelem,
                               optstruct *fpopts);

/** @cond */
#define FUNSUFFIX f
#define FPTYPE float
#define INTTYPE uint32_t
#define INTSUFFIX  U

#define DEFPREC   24
#define DEFEMAX  127
#define DEFEMIN -126
#define NLEADBITS  9
#define NBITS     32
#define FULLMASK 0xFFFFFFFFU
#define ABSMASK  0x7FFFFFFFU
#define SIGNMASK 0x80000000U
#define EXPMASK  0x7F800000U
#define FRACMASK 0x007FFFFFU

#ifdef PCG_VARIANTS_H_INCLUDED
#define MAXRAND 0xFFFFFFFFU
#define INITRAND(seed) pcg32_srandom_r(seed, time(NULL), (intptr_t)seed);
#define ADVANCERAND(seed, thread, nloc)                                        \
  pcg32_advance_r(seed, thread * nloc - 1);
#define GENRAND(seed) pcg32_random_r(seed)
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#warning "The default C random number generator is being used."
#warning "Please compile with -include <path-to-pcg_variants.h>"
#warning "and link with -L <path-to-libpcg_random.a> -lpcg_random."
#define MAXRAND 0x7FFFFFFFU
#define INITRAND(seed) *seed = time(NULL);
#define GEN_SINGLE_RAND(seed) ((INTTYPE)rand_r((unsigned int *)seed))
#endif /* #ifndef PCG_VARIANTS_H_INCLUDED */

#include "cpfloat_threshold_binary32.h"
#include "cpfloat_template.h"
/** @endcond */

#endif /* #ifndef _CPFLOAT_BINARY32_ */

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
