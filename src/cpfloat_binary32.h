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
doc_cpf_bivariate(sum, \f$ X_i = A_i + B_i \f$, 24, 127)
static inline int cpf_addf(float *X, const float *A, const float *B,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(difference, \f$ X_i = A_i - B_i \f$, 24, 127)
static inline int cpf_subf(float *X, const float *A, const float *B,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(product, \f$ X_i = A_i \times B_i \f$, 24, 127)
static inline int cpf_mulf(float *X, const float *A, const float *B,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(ratio, \f$ X_i = A_i / B_i \f$, 24, 127)
static inline int cpf_divf(float *X, const float *A, const float *B,
                           const size_t numelem, optstruct *fpopts);

/* Trigonometric functions. */
doc_cpf_univariate(trigonometric cosine, \f$ X_i = \cos(A_i) \f$, 24, 127)
static inline int cpf_cosf(float *X, const float *A,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(trigonometric sine, \f$ X_i = \sin(A_i) \f$, 24, 127)
static inline int cpf_sinf(float *X, const float *A,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(trigonometric tangent, \f$ X_i = \tan(A_i) \f$, 24, 127)
static inline int cpf_tanf(float *X, const float *A,
                           const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(inverse trigonometric cosine,
                   \f$ X_i = \mathrm{acos}(A_i) \f$, 24, 127)
static inline int cpf_acosf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse trigonometric sine,
                   \f$ X_i = \mathrm{asin}(A_i) \f$, 24, 127)
static inline int cpf_asinf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse trigonometric tangent,
                   \f$ X_i = \mathrm{atan}(A_i) \f$, 24, 127)
static inline int cpf_atanf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(2-argument arctangent,
                  \f$ X_i = \mathrm{atan} (B_i / A_i) \f$, 24, 127)
static inline int cpf_atan2f(float *X, const float *A, const float *B,
                             const size_t numelem, optstruct *fpopts);

/* Hyperbolic functions. */
doc_cpf_univariate(hyperbolic cosine, \f$ X_i = \mathrm{cosh}(A_i) \f$, 24, 127)
static inline int cpf_coshf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(hyperbolic sine, \f$ X_i = \mathrm{sinh}(A_i) \f$, 24, 127)
static inline int cpf_sinhf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(hyperbolic tangent , \f$ X_i = \mathrm{tanh}(A_i) \f$, 24, 127)
static inline int cpf_tanhf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(inverse hyperbolic cosine,
                   \f$ X_i = \mathrm{arcosh}(A_i) \f$, 24, 127)
static inline int cpf_acoshf(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse hyperbolic sine,
                   \f$ X_i = \mathrm{arsinh}(A_i) \f$, 24, 127)
static inline int cpf_asinhf(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(inverse hyperbolic tangent,
                   \f$ X_i = \mathrm{artanh}(A_i) \f$, 24, 127)
static inline int cpf_atanhf(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);

/* Exponentiation and logarithmic functions. */
doc_cpf_univariate(exponential, \f$ X_i = \exp(A_i) \f$, 24, 127)
static inline int cpf_expf(float *X, const float *A,
                           const size_t numelem, optstruct *fpopts);

doc_cpf_frexp(24, 127)
static inline int cpf_frexpf(float *X, int *exp, const float *A,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(2, 24, 127)
static inline int cpf_ldexpf(float *X, const float *A, const int *exp,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(natural logarithm, \f$ X_i = \log(A_i) \f$, 24, 127)
static inline int cpf_logf(float *X, const float *A,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-10 logarithm, \f$ X_i = \log_{10}(A_i) \f$, 24, 127)
static inline int cpf_log10f(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_modf(24, 127)
static inline int cpf_modff(float *X, float *intpart, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-2 exponential, \f$ X_i = 2^{A_i} \f$, 24, 127)
static inline int cpf_exp2f(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(exp(x) - 1, \f$ X_i = \exp(A_i) - 1 \f$, 24, 127)
static inline int cpf_expm1f(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_ilogb(24, 127)
static inline int cpf_ilogbf(int *exp, const float *A,
                             const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(natural logarithm of number shifted by one,
                   \f$ X_i = \log(1+A_i) \f$, 24, 127)
static inline int cpf_log1pf(float *X, const float *A,
                             size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-2 logarithm, \f$ X_i = \log_2(A_i) \f$, 24, 127)
static inline int cpf_log2f(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(base-FLT_RADIX logarithm of absolute value,
                   \f$ X_i = \log(\lvert A_i \rvert) \f$, 24, 127)
static inline int cpf_logbf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(FLT\_RADIX, 24, 127)
static inline int cpf_scalbnf(float *X, const float *A, const int *exp,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_scaling(FLT\_RADIX, 24, 127)
static inline int cpf_scalblnf(float *X, const float *A,
                               const long int *exp, const size_t numelem,
                               optstruct *fpopts);

/* Power functions. */
doc_cpf_bivariate(real powers, \f$ X_i = A_i^{B_i} \f$, 24, 127)
static inline int cpf_powf(float *X, const float *A, const float *B,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(square root, \f$ X_i = \sqrt{A_i} \f$, 24, 127)
static inline int cpf_sqrtf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(cube root, \f$ X_i = \sqrt[3]{A_i} \f$, 24, 127)
static inline int cpf_cbrtf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(hypotenuse of a right-angle triangle,
                  \f$ X_i = \sqrt{A_i^2 + B_i^2} \f$, 24, 127)
static inline int cpf_hypotf(float *X, const float *A, const float *B,
                             const size_t numelem, optstruct *fpopts);

/* Error and gamma functions. */
doc_cpf_univariate(error function, \f$ X_i = \mathrm{erf}(A_i) \f$, 24, 127)
static inline int cpf_erff(float *X, const float *A,
                           const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(complementary error function,
                   \f$ X_i = \mathrm{erfc}(A_i) \f$, 24, 127)
static inline int cpf_erfcf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(gamma function, \f$ X_i = \Gamma(A_i) \f$, 24, 127)
static inline int cpf_tgammaf(float *X, const float *A,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(natural logarithm of absolute value of gamma function,
                   \f$ X_i = \log(\lvert \Gamma(A_i) \rvert) \f$, 24, 127)
static inline int cpf_lgammaf(float *X, const float *A,
                              const size_t numelem, optstruct *fpopts);

/* Rounding and remainder functions. */
doc_cpf_univariate(ceiling function, \f$ X_i = \lceil A_i \rceil \f$, 24, 127)
static inline int cpf_ceilf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(floor function, \f$ X_i = \lfloor A_i \rfloor \f$, 24, 127)
static inline int cpf_floorf(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(floating-point remainder of division,
                  \f$ X_i = A_i \;\mathrm{mod}\; B_i \f$, 24, 127)
static inline int cpf_fmodf(float *X, const float *A, const float *B,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(integer truncation, \f$ X_i = \mathrm{trunc}(A_i) \f$, 24, 127)
static inline int cpf_truncf(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);

doc_cpf_univariate(closest integer (with round-to-nearest),
                   \f$ X_i = \mathrm{round}(A_i) \f$, 24, 127)
static inline int cpf_roundf(float *X, const float *A,
                             const size_t numelem, optstruct *fpopts);
doc_cpf_univariate(closest integer (with round-to-nearest),
                   \f$ X_i = \mathrm{round}(A_i) \f$, 24, 127)
static inline int cpf_lroundf(long *X, const float *A,
                              const size_t numelem, optstruct *fpopts);
doc_cpf_univariate_nobitflip(closest integer (with round-to-nearest),
                             \f$ X_i = \mathrm{round}(A_i) \f$, 24, 127)
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
                  \f$ X_i = A_i^2 - k \times B_i \f$
                  for largest \f$ k \f$ such that \f$ k \times B_i < A_i \f$,
                  24, 127)
static inline int cpf_remainderf(float *X, const float *A, const float *B,
                                 const size_t numelem, optstruct *fpopts);

doc_cpf_remquo(PMAX, EMAX)
static inline int cpf_remquof(float *X, int *quot,
                              const float *A, const float *B,
                              const size_t numelem, optstruct *fpopts);

/* Floating-point manipulation functions. */
doc_cpf_bivariate(number from magnitude and sign,
                  \f$ X_i = \mathrm{sign}(A_i) \times \lvert B_i \rvert \f$,
                  24, 127)
static inline int cpf_copysignf(float *X, const float *A, const float *B,
                                const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(next floating-point number in specified direction,
                  the floating-point number closest to \f$ A_i \f$ in the
                  direction of \f$ B_i \f$, 24, 127)
static inline int cpf_nextafterf(float *X, const float *A, const float *B,
                                 const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(next floating-point number in specified direction,
                  the floating-point number closest to \f$ A_i \f$ in the
                  direction of \f$ B_i \f$, 24, 127)
static inline int cpf_nexttowardf(float *X, const float *A,
                                  const long double *B,
                                  const size_t numelem,
                                  optstruct *fpopts);

/* Minimum, maximum, difference functions. */
doc_cpf_bivariate(positive difference, \f$ X_i = \lvert A_i - B_i \rvert \f$,
                  24, 127)
static inline int cpf_fdimf(float *X, const float *A, const float *B,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(element-wise maximum, \f$ X_i = \mathrm{max}(A_i, B_i) \f$,
                  24, 127)
static inline int cpf_fmaxf(float *X, const float *A, const float *B,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_bivariate(element-wise minimum, \f$ X_i = \mathrm{min}(A_i, B_i) \f$,
                  24, 127)
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
doc_cpf_univariate(absolute value, \f$ X_i = \lvert A_i \rvert \f$, 24, 127)
static inline int cpf_fabsf(float *X, const float *A,
                            const size_t numelem, optstruct *fpopts);
doc_cpf_trivariate(fused multiply-add , \f$ X_i = A_i \times B_i + C_i \f$,
                   24, 127)
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
