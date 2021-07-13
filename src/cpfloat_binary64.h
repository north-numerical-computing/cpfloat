/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/**
 * @file cpfloat_binary64.h
 * @brief CPFloat functions for rounding `double` arrays.
 *
 * @details This header file defines the function cpfloat() for rounding
 * `double` arrays to lower precision. When OpenMP support is specified at
 * compile time, the two functions cpfloat_sequential() and
 * cpfloat_parallel() are also defined.
 */
#ifndef _CPFLOAT_BINARY64_
#define _CPFLOAT_BINARY64_

#include "cpfloat_definitions.h"

/** @cond */
#define FUNSUFFIX
#define FPTYPE double
#define INTTYPE uint64_t
#define INTSUFFIX  ULL
#define DEFPREC    53
#define DEFEMAX  1023
#define DEFEMIN -1022
#define NLEADBITS  12
#define NBITS      64
#define FULLMASK 0xFFFFFFFFFFFFFFFFULL
#define ABSMASK  0x7FFFFFFFFFFFFFFFULL
#define SIGNMASK 0x8000000000000000ULL
#define EXPMASK  0x7FF0000000000000ULL
#define FRACMASK 0x000FFFFFFFFFFFFFULL

#ifdef PCG_VARIANTS_H_INCLUDED
#define NRNDBITS 64
#define SEEDTYPE pcg64_random_t
#define INITRAND_SINGLE(seed)                   \
  pcg64_srandom_r(seed,                         \
                  time(NULL),                   \
                  (intptr_t)seed)
#define INITRAND_MULTI(seed)                                    \
  pcg64_srandom_r(seed,                                         \
                  omp_get_thread_num() * 12345 + time(NULL),    \
                  (intptr_t)seed)
#define GENRAND(seed) pcg64_random_r(seed)
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#warning "The default C random number generator is being used."
#warning "Please compile with the option --include <path-to-pcg_variants.h>."
#define NRNDBITS 62
#define SEEDTYPE size_t
#define INITRAND_SINGLE(seed) *seed = time(NULL)
#define INITRAND_MULTI(seed) *seed = omp_get_thread_num() * 13254 + time(NULL)
#define GENRAND(seed) ((INTTYPE)rand_r((unsigned int *)seed) + \
                       ((INTTYPE)rand_r((unsigned int *)seed) << 31))
#endif /* #ifdef PCG_VARIANTS_H_INCLUDED */
/** @endcond */

/**
 * @brief Validate fields of @ref optstruct struct for `double` storage format.
 *
 * @details This function checks whether the parameters stored in @p fpopts are
 * valid when `double` is used as storage format.
 *
 * @param[in] fpopts Parameters describing the target format, the rounding mode,
 * and the probability of soft errors striking the rounded values.
 *
 * @return The function returns @b 0 if all the parameters are valid, and a
 * positive number if at least one of them is not. A negative number should be
 * understood as a warning, and indicates that cpfloat(),
 * cpfloat_sequential(), and cpfloat_parallel() will return @p 0 if @p fpopts
 * is used as fourth argument, but might not perform as intended.
 *
 * Possible return values are:
 * + @b -4 The rounding mode specified in @p fpopts->round does not correspond
 * to a valid choice, thus no rounding will be performed.
 * + @b -2 The required number of digits in @p fpopts->precision is between 26
 * and 53 inclusive, which might cause double rounding if round-to-nearest is
 * used.
 * + @b -1 The string in @p fpopts->format is not valid. This is not an error as
 * this value is not used by the C functions, but only by the MEX interface.
 * + @b  0 All the parameters in @p fpopts are valid.
 * + @b  2 The required number of digits in @p fpopts->precision is larger than
 * 53, the number of significant digits in a variable of type `double`.
 * + @b  3 The required maximum exponent in @p fpopts->emax is larger than 1023,
 * the largest possible exponent for a variable of type `double`.
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
int cpfloat_validate_optstruct(const optstruct *fpopts);

/**
 * @brief Round `double` array to lower precision.
 *
 * @details If the function executes without errors, then the array @p X
 * contains the @p numelem entries of the array @p A rounded to a
 * lower-precision target format. The parameters of the target format and the
 * rounding mode to be used are encoded in @p fpopts. If required, the function
 * flips one bit in some of the entries of @p X.
 *
 * If OpenMP support is specified when compiling the file that includes
 * cpfloat_binary64.h, then the function calls cpfloat_sequential() if @p
 * numelem is below the threshold parameter #OPENMP_THRESHOLD_double and
 * cpfloat_parallel() otherwise. When OpenMP is not supported, these two
 * functions are not defined and cpfloat() implements only the sequential
 * algorithm.
 *
 * @param[out] X Array of rounded values.
 * @param[in] A Input array.
 * @param[in] numelem Number of elements in @p X and @p A.
 * @param[in] fpopts Parameters describing the target format, the rounding mode,
 * and the probability of soft errors striking the rounded values.
 *
 * @return The function returns @b 1 if @p fpopts->precision is larger than 53,
 * @b 2 if @p fpopts->emax is larger than 1023, and @b 0 otherwise.
 */
static inline
int cpfloat(double *X,
            const double *A,
            const size_t numelem,
            const optstruct *fpopts);


#ifdef _OPENMP
/**
 * @brief Round `double` array to lower precision without using OpenMP.
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
 * @return The function returns @b 1 if @p fpopts->precision is larger than 53,
 * @b 2 if @p fpopts->emax is larger than 1023, and @b 0 otherwise.
 */
static inline
int cpfloat_sequential(double *X,
                       const double *A,
                       const size_t numelem,
                       const optstruct *fpopts);
/**
 * @brief Round `double` array to lower precision using multiple OpenMP threads.
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
 * @return The function returns @b 1 if @p fpopts->precision is larger than 53,
 * @b 2 if @p fpopts->emax is larger than 1023, and @b 0 otherwise.
 */
static inline
int cpfloat_parallel(double *X,
                     const double *A,
                     const size_t numelem,
                     const optstruct *fpopts);
#include "cpfloat_threshold_binary64.h"
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

#endif // _CPFLOAT_BINARY64_

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
