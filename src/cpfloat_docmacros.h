/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

#ifndef _CPFLOAT_DOCMACROS_
#define _CPFLOAT_DOCMACROS_

#define doc_cpfloat_validate_optstruct(FPTYPE, PMIN, PMAX, EMAX) \
/** \
 @brief Validate fields of @ref optstruct struct for `FPTYPE` storage format. \
 \
 @details This function checks whether the parameters stored in @p fpopts are \
 valid when `FPTYPE` is used as storage format. \
 \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 0 if all the parameters are valid, and a \
 positive number if at least one of them is not. A negative number should be \
 understood as a warning, and indicates that a CPFloat function @em will \
 return @p 0 if @p fpopts is used as fourth argument, but might not perform as \
 intended.<p/>\
 \
 Possible return values are: \
 \li @b -4 The rounding mode specified in @p fpopts->round does not correspond \
 to a valid choice, thus no rounding will be performed. \
 \li @b -2 The required number of digits in @p fpopts->precision is between \
 PMIN and PMAX inclusive, which might cause double rounding if round-to-\
 nearest is used. \
 \li @b -1 The string in @p fpopts->format is not valid. This is not an error \
 as this value is not used by the C functions, but only by the MEX interface. \
 \li @b  0 All the parameters in @p fpopts are valid. \
 \li @b  2 The required number of digits in @p fpopts->precision is larger \
 than PMAX, the number of significant digits in a variable of type `FPTYPE`. \
 \li @b  3 The required maximum exponent in @p fpopts->emax is larger than \
 EMAX, the largest possible exponent for a variable of type `FPTYPE`. \
 \li @b  5 The value of @p fpopts->flip indicates that soft errors should be \
 introduced, but @p fpopts->p is not a real number between 0 and 1 and thus \
 does not represent a valid probability.<p/>\
 \
 Errors take precedence over warnings, thus a nonpositive return value \
 implies no errors. In case of multiple issues, the return value is that of \
 the first error (or warning, if no error is present) encountered in the order \
 given in the list above. \
 */

#define doc_cpfloat(FPTYPE, PMAX, EMAX) \
/** \
 @brief Round `FPTYPE` array to lower precision. \
 \
 @details If the function executes without errors, then the array @p X \
 contains the @p numelem entries of the array @p A rounded to a \
 lower-precision target format. The parameters of the target format and the \
 rounding mode to be used are encoded in @p fpopts. If required, the function \
 flips one bit in some of the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine-dependent. <p/>\
 \
 @param[out] X Array of rounded values. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
 */

#define doc_cpf_univariate(MATHFUN, FUNSTRING, PMAX, EMAX) \
/** \
 @brief Compute MATHFUN rounded to lower precision. \
 \
 @details If the function executes without errors, then
 FUNSTRING \
 rounded to a lower-precision target format. The parameters of the \
 target format and the rounding mode to be used are encoded in @p fpopts. If \
 required, the function flips one bit in some of the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of rounded values. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_univariate_nobitflip(MATHFUN, FUNSTRING, PMAX, EMAX) \
/** \
 @brief Compute MATHFUN in lower precision. \
 \
 @details If the function executes without errors, then
 FUNSTRING \
 rounded to a lower-precision target format. The parameters of the \
 target format and the rounding mode to be used are encoded in @p fpopts.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of rounded values. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_bivariate(MATHFUN, FUNSTRING, PMAX, EMAX) \
/** \
 @brief Compute MATHFUN in lower precision. \
 \
 @details If the function executes without errors, then \
  FUNSTRING \
 rounded to a lower-precision target format. The parameters of the \
 target format and the rounding mode to be used are encoded in @p fpopts. If \
 required, the function flips one bit in some of the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of rounded values. \
 @param[in] A Input array. \
 @param[in] B Input array. \
 @param[in] numelem Number of elements in @p X, @p A, and @p B. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_trivariate(MATHFUN, FUNSTRING, PMAX, EMAX) \
/** \
 @brief Compute MATHFUN in lower precision. \
 \
 @details If the function executes without errors, then \
 FUNSTRING \
 rounded to a lower-precision target format. The parameters of the \
 target format and the rounding mode to be used are encoded in @p fpopts. If \
 required, the function flips one bit in some of the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of rounded values. \
 @param[in] A Input array. \
 @param[in] B Input array. \
 @param[in] C Input array. \
 @param[in] numelem Number of elements in @p X, @p A, @p B, and @p C. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_frexp(PMAX, EMAX) \
/** \
 @brief Exponent and normalized fraction of rounded floating-point number. \
 \
 @details If the function executes without errors, then: \
 \li if \f$ A_i \f$ is 0, then \f$ X_i \f$ and \f$ \exp_i \f$ are both set to \
 zero;\
 \li otherwise, \f$ X_i \f$ is a value in the range \f$ (-1;-0.5] \cup
 [0.5; 1) \f$ and \f$ \exp_i \f$ is an integer such that \f$ 2^{\exp_i} \
 \times X_i \f$ is equal to \f$ A_i \f$ rounded to a lower-precision target \
 format.<p/>\
 \
 The parameters of the target format and the rounding mode to be used are \
 encoded in @p fpopts. If \ required, the function flips one bit in some of \
 the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of floating-point values in \
 \f$ (-1;-0.5] \f$, \f$ [0.5; 1) \f$. \
 @param[out] exp Array of integer exponents. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X, @p A, @p B, and @p C. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_scaling(BASE, PMAX, EMAX) \
/** \
 @brief Scale number by power of BASE in lower precision. \
 \
 @details If the function executes without errors, then \f$ X_i = A_i \times \
 \mathrm{BASE}^{\exp_i} \f$ rounded to a lower-precision target format. \
 The parameters of the target format and the rounding mode to be used are \
 encoded in @p fpopts. If required, the function flips one bit in some of the \
 entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of rounded values. \
 @param[in] A Input array. \
 @param[in] exp Array of integer exponents. \
 @param[in] numelem Number of elements in @p X, @p A, @p B, and @p C. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_modf(PMAX, EMAX) \
/** \
 @brief Compute integral and fractional part. \
 \
 @details If the function executes without errors, then \f$ X_i \f$ is a value \
 the range \f$ (-1,1) \f$ and \f$ \mathrm{intpart}_i \f$ is an integer such \
 that \f$ X_i + \mathrm{intpart}_i \f$ is equal to \f$ A_i \f$ rounded to a \
 lower-precision target format. The parameters of the target format and the \
 rounding mode to be used are encoded in @p fpopts. If required, the function \
 flips one bit in some of the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of floating-point values in (-1, 1). \
 @param[out] intpart Array of integer parts. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X, @p A, @p B, and @p C. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_ilogb(PMAX, EMAX) \
/** \
 @brief Compute integral part of the logarithm of the absolute value. \
 \
 @details If the function executes without errors, the integer \f$ \exp_i \f$ \
 is the exponent used internally to express the floating-point value \
 \f$ A_i \f$ rounded to a lower-precision target format. In other words, \
 \f$ X_i \f$ is equal to \f$ \mathrm{trunc}(\log_b^{\lvert A_i \rvert}) \f$ \
 where \f$ b = \mathrm{FLT\_RADIX} \f$ is typically 2. The parameters of the \
 target format and the rounding mode to be used are encoded in @p fpopts.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] exp Array of floating-point values in \f$ (-1, 1) \f$. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X, @p A, @p B, and @p C. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_rint(PMAX, EMAX) \
/** \
 @brief Compute the closest integer with specified rounding mode. \
 \
 @details If the function executes without errors, then \f$ X_i \f$ is the \
 integral part of \f$ A_i \f$ rounded to a lower-precision target format and \
 \f$ \mathrm{exception}_i \f$ is set to 0 if \f$ X_i \f$ is equal to \
 \f$ A_i \f$ and to FE_INEXACT otherwise. The parameters of the target format \
 and the rounding mode to be used are encoded in @p fpopts. If required, the \
 function flips one bit in some of the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of rounded values. \
 @param[out] exception Array of floating-point exceptions. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_nearbyint(PMAX, EMAX) \
/** \
 @brief Compute the closest integer with specified rounding mode. \
 \
 @details If the function executes without errors, then \f$ X_i \f$ is the \
 integral part of \f$ A_i \f$ rounded to lower-precision target format. \
 The parameters of the target format and the rounding mode to be used are \
 encoded in @p fpopts. If required, the function flips one bit in some of the \
 entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of rounded values. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/


#define doc_cpf_remquo(PMAX, EMAX) \
/** \
 @brief Compute reminder and quotient of rounded numbers. \
 \
 @details If the function executes without errors, then \f$ \mathrm{quot}_i \f$ \
 and \f$ X_i \f$ are the (integral) quotient and the reminder of the division \
 \f$ A_i / B_i \f$  with \f$ A_i \f$ and \f$ B_i \f$ rounded to a \
 lower-precision target format. The parameters of the target format and the \
 rounding mode to be used are encoded in @p fpopts. If required, the function \
 flips one bit in some of the entries of @p X.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] X Array of reminders. \
 @param[out] quot Array of quotients. \
 @param[in] A Input array. \
 @param[in] B Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_fpclassify(PMAX, EMAX) \
/** \
 @brief Categorize floating-point values. \
 \
 @details If the function executes without errors, then \f$ r_i \f$ has value: \
 \li FP_INFINITE, if \f$ A_i \f$ is finite in the lower-precising target format; \
 \li FP_NAN, if \f$ A_i \f$ is a NaN in the lower-precising target format; \
 \li FP_NORMAL, if \f$ A_i \f$ is normal in the lower-precising target format; \
 \li FP_SUBNORMAL, if \f$ A_i \f$ is subnormal in the lower-precising target format; and \
 \li FP_ZERO, if \f$ A_i \f$ is zero in the lower-precising target format. <p/> \
 The parameters of the target format and the rounding mode to be used are \
 encoded in @p fpopts.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] r Array of classes. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#define doc_cpf_isfun(STRING, PMAX, EMAX) \
/** \
 @brief Check whether value is STRING in lower precision target format. \
 \
 @details If the function executes without errors, then \f$ r_i \f$ is a \
 nonzero integral value if \f$ A_i \f$ is STRING in a lower-precision target \
 format, and zero otherwise. The parameters of the target format and the \
 rounding mode to be used are encoded in @p fpopts.<p/>\
 \
 If OpenMP support is specified at compile time, several OpenMP threads are \
 used if @p numelem is large enough. This parameter is machine dependent.\
 \
 @param[out] r Array of Boolean values. \
 @param[in] A Input array. \
 @param[in] numelem Number of elements in @p X and @p A. \
 @param[in] fpopts Parameters describing the target format, the rounding mode, \
 and the probability of soft errors striking the rounded values. \
 \
 @return The function returns @b 1 if @p fpopts->precision is larger than \
 PMAX, @b 2 if @p fpopts->emax is larger than EMAX, and @b 0 otherwise. \
*/

#endif  /* #ifndef _CPFLOAT_DOCMACROS_ */

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
