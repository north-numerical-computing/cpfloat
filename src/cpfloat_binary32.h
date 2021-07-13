/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/**
 * @file cpfloat_binary32.h
 * @brief CPFloat functions for rounding `float` arrays.
 *
 * @details This header file defines the function cpfloatf() for rounding
 * `float` arrays to lower precision. When OpenMP support is specified at
 * compile time, the two functions cpfloatf_sequential() and
 * cpfloatf_parallel() are also defined.
 */

#ifndef _CPFLOAT_BINARY32_
#define _CPFLOAT_BINARY32_

#include "cpfloat_definitions.h"

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
#define NRNDBITS 32
#define SEEDTYPE pcg32_random_t
#define INITRAND_SINGLE(seed)               \
  pcg32_srandom_r(seed,                     \
                  time(NULL),               \
                  (intptr_t)seed)
#define INITRAND_MULTI(seed)                                    \
  pcg32_srandom_r(seed,                                         \
                  omp_get_thread_num() * 13254 + time(NULL),    \
                  (intptr_t)seed)
#define GENRAND(seed) pcg32_random_r(seed)
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#warning "The default C random number generator is being used."
#warning "Please compile with the option --include <path-to-pcg_variants.h>."
#define NRNDBITS 31
#define SEEDTYPE size_t
#define INITRAND_SINGLE(seed) *seed = time(NULL)
#define INITRAND_MULTI(seed) *seed = omp_get_thread_num() * 13254 + time(NULL)
#define GENRAND(seed) (INTTYPE)rand_r((unsigned int *)seed)
#endif /* #ifdef PCG_VARIANTS_H_INCLUDED */
/** @endcond */

/**
 * @brief Validate fields of @ref optstruct struct for `float` storage format.
 *
 * @details This function checks whether the parameters stored in @p fpopts are
 * valid when `float` is used as storage format.
 *
 * @param[in] fpopts Parameters describing the target format, the rounding mode,
 * and the probability of soft errors striking the rounded values.
 *
 * @return The function returns @b 0 if all the parameters are valid, and a
 * positive number if at least one of them is not. A negative number should be
 * understood as a warning, and indicates that cpfloatf(),
 * cpfloatf_sequential(), and cpfloatf_parallel() will return @p 0 if @p fpopts
 * is used as fourth argument, but might not perform as intended.
 *
 * Possible return values are:
 * + @b -4 The rounding mode specified in @p fpopts->round does not correspond
 * to a valid choice, thus no rounding will be performed.
 * + @b -2 The required number of digits in @p fpopts->precision is between 12
 * and 24 inclusive, which might cause double rounding if round-to-nearest is
 * used.
 * + @b -1 The string in @p fpopts->format is not valid. This is not an error as
 * this value is not used by the C functions, but only by the MEX interface.
 * + @b  0 All the parameters in @p fpopts are valid.
 * + @b  2 The required number of digits in @p fpopts->precision is larger than
 * 24, the number of significant digits in a variable of type `float`.
 * + @b  3 The required maximum exponent in @p fpopts->emax is larger than 127,
 * the largest possible exponent for a variable of type `float`.
 * + @b  5 The value of @p fpopts->flip indicates that soft errors should be
 * introduced, but @p fpopts->p is not a real number between 0 and 1 and thus
 * does not represent a valid probability.
 *
 * Errors take precedence over warnings, thus a nonpositive return values
 * implies no errors. In case of multiple issues, the return value is that of
 * the first error (or warning, if no error is present) encountered in the order
 * given in the list above.
 */
static inline
int cpfloatf_validate_optstruct(const optstruct *fpopts);

/**
 * @brief Round `float` array to lower precision.
 *
 * @details If the function executes without errors, then the array @p X
 * contains the @p numelem entries of the array @p A rounded to a
 * lower-precision target format. The parameters of the target format and the
 * rounding mode to be used are encoded in @p fpopts. If required, the function
 * flips one bit in some of the entries of @p X.
 *
 * If OpenMP support is specified when compiling the file that includes
 * cpfloat_binary32.h, then the function calls cpfloatf_sequential() if @p
 * numelem is below the threshold parameter #OPENMP_THRESHOLD_float and
 * cpfloatf_parallel() otherwise. When OpenMP is not supported, these two
 * functions are not defined and cpfloatf() implements only the sequential
 * algorithm.
 *
 * @param[out] X Array of rounded values.
 * @param[in] A Input array.
 * @param[in] numelem Number of elements in @p X and @p A.
 * @param[in] fpopts Parameters describing the target format, the rounding mode,
 * and the probability of soft errors striking the rounded values.
 *
 * @return The function returns @b 1 if @p fpopts->precision is larger than 24,
 * @b 2 if @p fpopts->emax is larger than 127, and @b 0 otherwise.
 */
static inline
int cpfloatf(float *X,
             const float *A,
             const size_t numelem,
             const optstruct *fpopts);

#ifdef _OPENMP
/**
 * @brief Round `float` array to lower precision without using OpenMP.
 *
 * @details If the function executes without errors, then the array @p X
 * contains the @p numelem entries of the array @p A rounded to a
 * lower-precision target format. The parameters of the target format and the
 * rounding mode to be used are encoded in @p fpopts. If required, the function
 * flips one bit in some of the entries of @p X.
 *
 * @param[out] X Array of rounded values.
 * @param[in] A Input array.
 * @param[in] numelem Number of elements in @p X and @p A.
 * @param[in] fpopts Parameters for target format, rounding mode, and soft errors.
 *
 * @return The function returns @b 1 if @p fpopts->precision is larger than 24,
 * @b 2 if @p fpopts->emax is larger than 127, and @b 0 otherwise.
 */
static inline
int cpfloatf_sequential(float *X,
                        const float *A,
                        const size_t numelem,
                        const optstruct *fpopts);
/**
 * @brief Round `float` array to lower precision using multiple OpenMP threads.
 *
 * @details If the function executes without errors, then the array @p X
 * contains the @p numelem entries of the array @p A rounded to a
 * lower-precision target format. The parameters of the target format and the
 * rounding mode to be used are encoded in @p fpopts. If required, the function
 * flips one bit in some of the entries of @p X.
 *
 * This function tries to use as many OpenMP threads as available on the system.
 *
 * @param[out] X Array of rounded values.
 * @param[in] A Input array.
 * @param[in] numelem Number of elements in @p X and @p A.
 * @param[in] fpopts Parameters for target format, rounding mode, and soft errors.
 *
 * @return The function returns @b 1 if @p fpopts->precision is larger than 24,
 * @b 2 if @p fpopts->emax is larger than 127, and @b 0 otherwise.
 */
static inline
int cpfloatf_parallel(float *X,
                      const float *A,
                      const size_t numelem,
                      const optstruct *fpopts);
#include "cpfloat_threshold_binary32.h"
/** @cond */
#define USE_OPENMP
#include "cpfloat_template.h"
#undef USE_OPENMP
/** @endcond */
#endif /* #ifdef _OPENMP */

#define SINGLE_THREADED
/** @cond */
#include "cpfloat_template.h"
/** @endcond */
#undef SINGLE_THREADED

#endif // _CPFLOAT_BINARY32_

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
