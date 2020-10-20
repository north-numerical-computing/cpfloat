/**
 * @file cpfloat_threshold_binary64.h
 * @brief Size of smallest `double` array on which to use multiple OpenMP threads.
 */

/**
 * @brief Size of smallest array on which cpfloat() uses multiple threads.
 *
 * @details Threshold for switching between cpfloat_sequential() and
 * cpfloat_parallel() in cpfloat(). The value of this constant is ignored
 * if the file that includes cpfloat_binary64.h is compiled without OpenMP
 * support.
 */
#define OPENMP_THRESHOLD_double 1
