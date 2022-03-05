/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/**
 * @file cpfloat_threshold_binary32.h
 * @brief Size of smallest `float` array on which to use multiple OpenMP threads.
 */

/**
 * @brief Size of smallest array on which cpfloatf() uses multiple threads.
 *
 * @details Threshold for switching between cpfloatf_sequential() and
 * cpfloatf_parallel() in cpfloatf(). The value of this constant is ignored
 * if the file that includes cpfloat_binary32.h is compiled without OpenMP
 * support.
 */
#define OPENMP_THRESHOLD_float 1
