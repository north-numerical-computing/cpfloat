/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * This file is part of CPFloat.
 *
 * The macros and functions defined here require a number of macros be defined
 * before this header is included. This file should never be included directly
 * in production code. In order to use CPFloat, it suffices to include
 * cpfloat_binary32.h and cpfloat_binary64.h. These two headers include this
 * file and handle the definition of macros correctly.
 */

/******************************************
 ******************************************
 * TYPE-INDEPENDENT MACROS AND FUNCTIONS. *
 ******************************************
 ******************************************/

#ifndef TYPE_INDEPENDENT_MACROS_INCLUDED
#define TYPE_INDEPENDENT_MACROS_INCLUDED

/* Define empty macros. */
#define NOOP
#define NOARG
#define UNUSED(x) (void)(x);

/* Remove one layer of parentheses, if any. */
#define DO_REMOVE_PARENTHESES
#define EXPAND_INTERIOR(...) REMOVE_INTERIOR_(__VA_ARGS__)
#define REMOVE_INTERIOR_(...) DO_ ## __VA_ARGS__
#define REMOVE_PARENTHESES(...) REMOVE_PARENTHESES __VA_ARGS__
#define DEPARENTHESIZE_MAYBE(X) EXPAND_INTERIOR(REMOVE_PARENTHESES X)

/* Concatenate strings in macros. */
#define CONCATENATE_INNER(arg1, arg2) arg1 ## arg2
#define CONCATENATE(arg1, arg2) CONCATENATE_INNER(arg1, arg2)

/* Functions to initialize bit pseudo-random generators and generate bits. */
#define BITSEED bitseed
#define BITSEEDTYPE cpfloat_bitseed_t
#define BITTYPE unsigned int

#ifdef PCG_VARIANTS_H_INCLUDED
#define INITBIT(seed) pcg32_srandom_r(seed, time(NULL), (intptr_t)seed);
#define ADVANCEBIT(seed, thread, nloc) \
  pcg32_advance_r(seed, thread * nloc - 1);
#define GENBIT(seed) (pcg32_random_r(seed) & (1U << 31))
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#define INITBIT(seed) *seed = time(NULL);
#define GEN_SINGLE_BIT(seed) (rand_r(seed) & (1U << 30))
#ifdef _OPENMP
#define PRNG_ADVANCE_BIT prng_advance_bit
static inline BITTYPE PRNG_ADVANCE_BIT(BITSEEDTYPE *seed, size_t delta) {
  for (size_t i=0; i<delta; i++)
    rand_r(seed);
  return GEN_SINGLE_BIT(seed);
}
#define ADVANCEBIT(seed, thread, nloc) PRNG_ADVANCE_BIT(seed, thread);
#define GENBIT(seed) PRNG_ADVANCE_BIT(seed, nthreads)
#else /* #ifdef _OPENMP */
#define GENBIT(seed) (GEN_SINGLE_BIT(seed))
#endif /* #ifdef _OPENMP */
#endif  /* #ifdef PCG_VARIANTS_H_INCLUDED */

#define PRNG_BIT_INIT                                                          \
  if (fpopts->BITSEED == NULL) {                                               \
    fpopts->BITSEED= malloc(sizeof(*fpopts->BITSEED));                         \
    INITBIT(fpopts->BITSEED)                                                   \
  }

#define PRNG_BIT_ADVANCE(PARALLEL) PRNG_BIT_ADVANCE_##PARALLEL
#define PRNG_BIT_ADVANCE_SEQ                                                   \
  size_t nthreads = 1;                                                         \
  UNUSED(nthreads)                                                             \
  size_t istart = 0;                                                           \
  size_t iend = numelem;                                                       \
  params.BITSEED = fpopts->BITSEED;
#define PRNG_BIT_ADVANCE_PAR                                                   \
  size_t nthreads = omp_get_num_threads();                                     \
  size_t nloc = ceil((double)numelem / nthreads);                              \
  size_t thread = omp_get_thread_num();                                        \
  size_t istart = nloc * thread;                                               \
  size_t iend = nloc * (thread + 1);                                           \
  iend = iend > numelem ? numelem : iend;                                      \
  BITSEEDTYPE tmp_bitseed = *(fpopts->BITSEED);                                \
  params.BITSEED = &tmp_bitseed;                                               \
  ADVANCEBIT(params.BITSEED, thread, nloc)

#define PRNG_BIT_UPDATE(PARALLEL) PRNG_BIT_UPDATE_##PARALLEL
#define PRNG_BIT_UPDATE_SEQ
#define PRNG_BIT_UPDATE_PAR                                                    \
  if (thread == nthreads - 1) {                                                \
    *(fpopts->BITSEED) = *(params.BITSEED);                                    \
  }

/*
 * Define type-independent functions.
 * Prototypes are in cpfloat_definitions.
 */

optstruct *init_optstruct() {
  optstruct *fpopts = malloc(sizeof(*fpopts));
  fpopts->bitseed = NULL;
  fpopts->randseedf = NULL;
  fpopts->randseed = NULL;
  return fpopts;
}

int free_optstruct(optstruct *fpopts) {
  if (fpopts == NULL)
    return -1;
  else {
    if (fpopts->bitseed != NULL)
      free(fpopts->bitseed);
    if (fpopts->randseedf != NULL)
      free(fpopts->randseedf);
    if (fpopts->randseed != NULL)
      free(fpopts->randseed);
    free(fpopts);
    return 0;
  }
}

#endif /* #ifndef TYPE_INDEPENDENT_MACROS_INCLUDED */



/****************************************
 ****************************************
 * TYPE-DEPENDENT MACROS AND FUNCTIONS. *
 ****************************************
 ****************************************/

#define ADDSUFFIXTO(x) CONCATENATE(x, FUNSUFFIX)

/* Struct to reinterpret floating-point values as unsigned integers. */
#define FPUNION ADDSUFFIXTO(fpint)
typedef union {
  INTTYPE intval;
  FPTYPE fpval;
} FPUNION;

#define INTCONST(x) CONCATENATE(x, INTSUFFIX)

#define INTOF(x)(((FPUNION *)(x))->intval)
#define INTOFCONST(x)(((FPUNION){.fpval = (FPTYPE)x}).intval)
#define FPOF(x)(((FPUNION){.intval = (INTTYPE)(x)}).fpval)


/* Types for internal state of pseudo-random number generator. */
#define RANDSEED ADDSUFFIXTO(randseed)
#define RANDSEEDTYPE CONCATENATE(ADDSUFFIXTO(cpfloat_randseed),_t)

/* Relevant parameters of floating-point number system. */
#define FPPARAMS ADDSUFFIXTO(fpparams)
typedef struct {
  cpfloat_precision_t precision;
  cpfloat_exponent_t emax;
  cpfloat_exponent_t emin;
  cpfloat_subnormal_t subnormal;
  cpfloat_rounding_t round;
  FPTYPE ftzthreshold;
  FPTYPE xmin;
  FPTYPE xmax;
  FPTYPE xbnd;
  INTTYPE leadmask;
  INTTYPE trailmask;
  BITSEEDTYPE *BITSEED;
  RANDSEEDTYPE *RANDSEED;
} FPPARAMS;

/* Parameters used for subnormal numbers. */
#define LOCPARAMS ADDSUFFIXTO(locparams)
typedef struct {
  cpfloat_precision_t precision;
  INTTYPE leadmask;
  INTTYPE trailmask;
} LOCPARAMS;

#define SIGN(x)(SIGNMASK & INTOF(x))
#define ABS(x)(FPOF((INTTYPE)(ABSMASK & INTOF(x))))

#ifndef PCG_VARIANTS_H_INCLUDED
#ifdef _OPENMP
#define PRNG_ADVANCE_RAND ADDSUFFIXTO(prng_advance_rand)
static inline INTTYPE PRNG_ADVANCE_RAND(RANDSEEDTYPE *seed, size_t delta) {
  for (size_t i=0; i<delta; i++)
    rand_r((unsigned int *)seed);
  return GEN_SINGLE_RAND(seed);
}
#define ADVANCERAND(seed, thread, nloc) PRNG_ADVANCE_RAND(seed, thread);
#define GENRAND(seed) PRNG_ADVANCE_RAND(seed, nthreads)
#else /* #ifdef _OPENMP */
#define GENRAND(seed) (GEN_SINGLE_RAND(seed))
#endif /* #ifdef _OPENMP */
#endif  /* #ifdef PCG_VARIANTS_H_INCLUDED */

/* Functions to initialize the pseudo-random number generator. */
#define PRNG_RAND_INIT                                                         \
  if (fpopts->RANDSEED == NULL) {                                              \
    fpopts->RANDSEED = malloc(sizeof(*fpopts->RANDSEED));                      \
    INITRAND(fpopts->RANDSEED)                                                 \
  }

#define PRNG_RAND_ADVANCE(PARALLEL) PRNG_RAND_ADVANCE_##PARALLEL
#define PRNG_RAND_ADVANCE_SEQ                                                  \
  size_t nthreads = 1;                                                         \
  UNUSED(nthreads)                                                             \
  size_t istart = 0;                                                           \
  size_t iend = numelem;                                                       \
  params.RANDSEED = fpopts->RANDSEED;
#define PRNG_RAND_ADVANCE_PAR                                                  \
  size_t nthreads = omp_get_num_threads();                                     \
  size_t nloc = ceil((double)numelem / nthreads);                              \
  size_t thread = omp_get_thread_num();                                        \
  size_t istart = nloc * thread;                                               \
  size_t iend = nloc * (thread + 1);                                           \
  iend = iend > numelem ? numelem : iend;                                      \
  RANDSEEDTYPE tmp_randseed = *(fpopts->RANDSEED);                             \
  params.RANDSEED = &tmp_randseed;                                             \
  ADVANCERAND(params.RANDSEED, thread, nloc)

#define PRNG_RAND_UPDATE(PARALLEL) PRNG_RAND_UPDATE_##PARALLEL
#define PRNG_RAND_UPDATE_SEQ
#define PRNG_RAND_UPDATE_PAR                                                   \
  if (thread == nthreads - 1) {                                                \
    *(fpopts->RANDSEED) = *(params.RANDSEED);                                  \
  }

/* Function to validate floating-point options passed to CPFloat interfaces. */
#define VALIDATE_INPUT ADDSUFFIXTO(CONCATENATE(cpfloat,_validate_optstruct))
static inline int VALIDATE_INPUT(const optstruct *fpopts) {

  int retval;

  /* Set retval to -1 if format is not valid. */
  const char *valid_formats [] = {"",
                                  "bfloat", "bf16", "b",
                                  "binary16", "fp16", "half", "h",
                                  "TensorFloat-32", "tf32", "t",
                                  "binary32", "fp32", "single", "s",
                                  #if DEFPREC >= 53
                                  "binary64", "fp64", "double", "d",
                                  #endif /* #if DEFPREC >= 53 */
                                  "custom", "c"};
  #if DEFPREC >= 53
  size_t nformats = 21;
  #else /* #if DEFPREC >= 53 */
  size_t nformats = 17;
  #endif /* #if DEFPREC >= 53 */
  retval = -1;
  for (size_t i=0; i<nformats; i++) {
    if (strcmp(fpopts->format, valid_formats[i]))
      retval = 0;
  }

  /* Return 2 if precision is invalid (either nonpositive or too large). */
  if (fpopts->precision > DEFPREC || fpopts->precision <= 0)
    return 2;

  /* Set retval to -2 if double rounding may occur. */
  if (fpopts->precision > (fpopts->round<=1 ? floor((DEFPREC-2)/2) : DEFPREC-1))
    retval = -2;

  /* Return 3 if emax is invalid (either nonpositive or too large). */
  if (fpopts->emax > DEFEMAX || fpopts->emax <= 0)
    return 3;

  /* Set retval to -4 if rounding mode is set to no rounding. */
  if (fpopts->round == CPFLOAT_NO_RND)
    retval = -4;

  /* Return 5 if p is required but is not a valid probability. */
  if (fpopts->flip != CPFLOAT_NO_SOFTERR && (fpopts->p > 1 || fpopts->p < 0))
    return 5;

  /* Return 0 or warning value. */
  return retval;
}

/* Compute floating-point parameters required by the rounding functions. */
#define COMPUTE_GLOBAL_PARAMS ADDSUFFIXTO(compute_golbal_params)
static inline FPPARAMS COMPUTE_GLOBAL_PARAMS(const optstruct *fpopts,
                                             int *retval) {

  /* Actual precision and exponent range. */
  *retval = 0;
  cpfloat_precision_t precision = fpopts->precision;
  cpfloat_exponent_t emax = fpopts->explim == CPFLOAT_EXPRANGE_TARG ?
    fpopts->emax :
    DEFEMAX;
  if (precision > DEFPREC) {
    precision = DEFPREC;
    *retval = 1;
  }
  if (emax > DEFEMAX) {
    emax = DEFEMAX;
    *retval = 2;
  }

  /* Derived floating point parameters. */
  int emin = 1-emax;
  FPTYPE xmin = ldexp(1., emin);              /* Smallest pos. normal. */
  FPTYPE xmins = ldexp(1., emin-precision+1); /* Smallest pos. subnormal. */
  FPTYPE ftzthreshold = (fpopts->subnormal == CPFLOAT_SUBN_USE) ? xmins : xmin;
  FPTYPE xmax = ldexp(1., emax) * (2-ldexp(1., 1-precision));
  FPTYPE xbnd = ldexp(1., emax) * (2-ldexp(1., -precision));

  /* Bitmasks. */
  INTTYPE leadmask = FULLMASK << (DEFPREC-precision); /* To keep. */
  INTTYPE trailmask = leadmask ^ FULLMASK;            /* To discard. */

  FPPARAMS params = {precision, emax, emin, fpopts->subnormal, fpopts->round,
                     ftzthreshold, xmin, xmax, xbnd,
                     leadmask, trailmask, NULL, NULL};

  return params;
}

/* Compute floating point parameters required for rounding subnormals. */
#define UPDATE_GLOBAL_BITMASKS ADDSUFFIXTO(update_global_bitmasks)
static inline void UPDATE_GLOBAL_BITMASKS(FPPARAMS *params) {

  /* Bitmasks. */
  params->leadmask = FULLMASK << (DEFPREC-params->precision); /* To keep. */
  params->trailmask = params->leadmask ^ FULLMASK;            /* To discard. */
}

/* Compute floating point parameters required for rounding subnormals. */
#define UPDATE_LOCAL_PARAMS ADDSUFFIXTO(update_local_params)
static inline void UPDATE_LOCAL_PARAMS(const FPTYPE *A,
                                       const FPPARAMS *params,
                                       LOCPARAMS *lparams) {
  if (ABS(A) < params->xmin && params->emin > 1-DEFEMAX) {
    int ndigits = params->precision - params->emin - DEFEMAX +
      ((EXPMASK & INTOF(A)) >> (DEFPREC - 1));
    lparams->precision = ndigits > 0 ? ndigits : 0;
    INTTYPE leadmask = params->leadmask;
    leadmask ^= (((INTCONST(1) <<
                   (params->precision - lparams->precision)) - INTCONST(1))
                 << (DEFPREC-params->precision));
    lparams->leadmask = leadmask;
    lparams->trailmask = leadmask ^ FULLMASK;
  } else {
    lparams->precision = params->precision;
    lparams->leadmask = params->leadmask;
    lparams->trailmask = params->trailmask;
  }
}

/*******************************
 * Macros for scalar rounding. *
 *******************************/

/*
 * These macros round scalars of type FPTYPE. The input arguments are:
 *   + FPTYPE *x, address where to store the output;
 *   + FPTYPE *y, address where the input is stored;
 *   + FPPARAMS *p, address where the target format parameters are stored;
 *   + LOCPARAMS *lp, address where the parameters for subnormals are stored.
 *
 * For some of the rounding modes, it is possible to simplify the algorithms
 * when storage and target format have the same exponent range. For each of
 * these rounding modes, we provide two macros: one with the suffix _SAME_EXP,
 * which assumes that storage and target format dedicate the same number of bits
 * to the exponent field, and one with the suffix _OTHER_EXP, which only assumes
 * that the exponent range of the target format is a subrange of that of the
 * storage format.
 *
 */

/* Round-to-nearest with ties-to-away. */
#define RN_TIES_TO_AWAY_SCALAR_SAME_EXP(x, y, p)                               \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {                  \
    if (ABS(y) >= p->xmin/2)                                                   \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmin));                              \
    else                                                                       \
      *(x) = FPOF(SIGN(y));                                                    \
  } else                                                                       \
    *(x) = FPOF((INTOF(y) +                                                    \
                 (INTCONST(1) << (DEFPREC-1-p->precision))) & p->leadmask);

#define RN_TIES_TO_AWAY_SCALAR_OTHER_EXP(x, y, p, lp)                          \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    if (ABS(y) < p->ftzthreshold/2)                                            \
      *(x) = FPOF(SIGN(y));                                                    \
    else                                                                       \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));                      \
  } else if (ABS(y) >= p->xbnd) { /* Overflow */                               \
    *(x) = FPOF(SIGN(y) | INTOFCONST(INFINITY));                               \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    *(x) = FPOF((INTOF(y) +                                                    \
                 (INTCONST(1) << (DEFPREC-1-lp->precision))) & lp->leadmask);  \
  }

/* Round-to-nearest with ties-to-zero. */
#define RN_TIES_TO_ZERO_SCALAR_SAME_EXP(x, y, p)                               \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {                  \
    if (ABS(y) > p->xmin/2)                                                    \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmin));                              \
    else                                                                       \
      *(x) = FPOF(SIGN(y));                                                    \
  } else                                                                       \
    *(x) = FPOF((INTOF(y) + (p->trailmask>>1)) & p->leadmask);

#define RN_TIES_TO_ZERO_SCALAR_OTHER_EXP(x, y, p, lp)                          \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    if (ABS(y) > p->ftzthreshold/2)                                            \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));                      \
    else                                                                       \
      *(x) = FPOF(SIGN(y));                                                    \
  } else if (ABS(y) > p->xbnd) { /* Overflow */                                \
    *(x) = FPOF(SIGN(y) | INTOFCONST(INFINITY));                               \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    *(x) = FPOF((INTOF(y) + (lp->trailmask>>1)) & lp->leadmask);               \
  }

/* Round-to-nearest with ties-to-even. */
#define RN_TIES_TO_EVEN_SCALAR_SAME_EXP(x, y, p)                               \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {                  \
    if (ABS(y) <= p->xmin/2)                                                   \
      *(x) = FPOF(SIGN(y));                                                    \
    else                                                                       \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmin));                              \
  } else {                                                                     \
    INTTYPE LSB = (INTOF(y) >> (DEFPREC-p->precision)) & INTCONST(1);          \
    *(x) = FPOF((INTOF(y) + (p->trailmask >> 1) + LSB) & p->leadmask);         \
  }

#define RN_TIES_TO_EVEN_SCALAR_OTHER_EXP(x, y, p, lp)                          \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    if (ABS(y) <= p->ftzthreshold/2)                                           \
      *(x) = FPOF(SIGN(y));                                                    \
    else                                                                       \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));                      \
  } else if (ABS(y) >= p->xbnd) { /* Overflow */                               \
    *(x) = FPOF(SIGN(y) | INTOFCONST(INFINITY));                               \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    INTTYPE LSB = ((INTOF(y) >> (DEFPREC-lp->precision)) & INTCONST(1))        \
      | (lp->precision == 1 && DEFEMAX != p->emax); /* Hidden bit is one. */   \
    *(x) = FPOF((INTOF(y) + ((lp->trailmask >> 1) + LSB)) & lp->leadmask);     \
  }

/* Round-to-plus-infinity (also known as round-up). */
#define RD_TWD_PINF_SCALAR_SAME_EXP(x, y, p)                                   \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {                  \
    if(*(y) > 0)                                                               \
      *(x) = p->xmin;                                                          \
    else                                                                       \
      *(x) = FPOF(SIGN(y));                                                    \
  } else {                                                                     \
    if (*(y) > 0) /* Add ulp if x is positive finite. */                       \
      *(x) = FPOF((INTOF(y) + p->trailmask) & p->leadmask);                    \
    else                                                                       \
      *(x) = FPOF(INTOF(y) & p->leadmask);                                     \
  }

#define RD_TWD_PINF_SCALAR_OTHER_EXP(x, y, p, lp)                              \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    *(x) = *(y) > 0 ? p->ftzthreshold : 0;                                     \
  } else if (ABS(y) > p->xmax) { /* Overflow */                                \
    if (*(y) > p->xmax)                                                        \
      *(x) = INFINITY;                                                         \
    else if (*(y) < -p->xmax && *(y) != -INFINITY)                             \
      *(x) = -p->xmax;                                                         \
    else /* *(y) == -INFINITY */                                               \
      *(x) = -INFINITY;                                                        \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    if (SIGN(y) == 0) /* Add ulp if x is positive.  */                         \
      *(x) = FPOF((INTOF(y) + lp->trailmask) & lp->leadmask);                  \
    else                                                                       \
      *(x) = FPOF(INTOF(y) & lp->leadmask);                                    \
  }

/* Round-to-minus-infinity (also known as round-down). */
#define RD_TWD_NINF_SCALAR_SAME_EXP(x, y, p)                                   \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {                  \
    if(*(y) < 0)                                                               \
      *(x) = -p->xmin;                                                         \
    else                                                                       \
      *(x) = FPOF(SIGN(y));                                                    \
  } else {                                                                     \
    if (*(y) < 0) /* Subtract ulp if x is positive finite. */                  \
      *(x) = FPOF((INTOF(y) + p->trailmask) & p->leadmask);                    \
    else                                                                       \
      *(x) = FPOF(INTOF(y) & p->leadmask);                                     \
  }

#define RD_TWD_NINF_SCALAR_OTHER_EXP(x, y, p, lp)                              \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    *(x) = *(y) >= 0 ? 0 : -p->ftzthreshold;                                   \
  } else if (ABS(y) > p->xmax) { /* Overflow */                                \
    if (*(y) < -p->xmax)                                                       \
      *(x) = -INFINITY;                                                        \
    else if (*(y) > p->xmax && *(y) != INFINITY)                               \
      *(x) = p->xmax;                                                          \
    else /* *(y) == INFINITY */                                                \
      *(x) = INFINITY;                                                         \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    if (SIGN(y)) /* Subtract ulp if x is positive. */                          \
      *(x) = FPOF((INTOF(y) + lp->trailmask) & lp->leadmask);                  \
    else                                                                       \
      *(x) = FPOF(INTOF(y) & lp->leadmask);                                    \
  }

/* Round-to-zero (also known as truncation). */
#define RD_TWD_ZERO_SCALAR_SAME_EXP(x, y, p)                                   \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin)                    \
    *(x) = FPOF(SIGN(y));                                                      \
  else                                                                         \
    *(x) = FPOF(INTOF(y) & p->leadmask);

#define RD_TWD_ZERO_SCALAR_OTHER_EXP(x, y, p, lp)                              \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    *(x) = FPOF(SIGN(y));                                                      \
  } else if (ABS(y) > p->xmax && ABS(y) != INFINITY) { /* Overflow */          \
    *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmax));                                \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    *(x) = FPOF(INTOF(y) & lp->leadmask);                                      \
  }

/* Stochastic rounding with proportional probabilities. */
#define RS_PROP_SCALAR(x, y, p, lp)                                            \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    int expx = ((EXPMASK & INTOF(y)) >> (DEFPREC - 1)) - DEFEMAX;              \
    INTTYPE trailfrac = (INTOF(y) & (FULLMASK >> NLEADBITS));                  \
    if (expx == -DEFEMAX) { /* No implicit bit (x is subnormal). */            \
      expx += 1;                                                               \
    } else { /* Make implicit bit explicit (x is normal). */                   \
      trailfrac = (INTOF(y) & (FULLMASK >> NLEADBITS))                         \
        | (INTCONST(1) << (DEFPREC-1));                                        \
    }                                                                          \
    int localemin = p->subnormal == CPFLOAT_SUBN_USE ?                         \
      p->emin - (int)p->precision + 1: p->emin;                                \
    int expdiff = localemin - expx;                                            \
    INTTYPE rnd = GENRAND(p->RANDSEED) >> 1;                                   \
    /* Shift fraction of *(y) left or right as needed. */                      \
    if (expdiff <= NLEADBITS - 1)                                              \
      trailfrac <<= NLEADBITS - 1 - expdiff;                                   \
    else {                                                                     \
      if (expdiff - (NLEADBITS - 1) >= NBITS)                                  \
        trailfrac = 0;                                                         \
      else                                                                     \
        trailfrac >>= (expdiff - (NLEADBITS - 1));                             \
    }                                                                          \
    if (trailfrac > rnd) {                                                     \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));                      \
    } else {                                                                   \
      *(x) = FPOF(SIGN(y));                                                    \
      continue;                                                                \
    }                                                                          \
  } else if (ABS(y) < p->xbnd) { /* Rounding possibly required. */             \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    INTTYPE rndbuf = GENRAND(p->RANDSEED) & lp->trailmask;                     \
    if (ABS(y) > p->xmax) {                                                    \
      rndbuf = rndbuf >> 1;                                                    \
      lp->leadmask = (lp->leadmask >> 1) | SIGNMASK;                           \
    }                                                                          \
    *(x) = FPOF((INTOF(y) + rndbuf) & lp->leadmask);                           \
  } else {                                                                     \
    *(x) = *(y);                                                               \
  }                                                                            \
  if (ABS(x) >= p->xbnd) /* Overflow */                                        \
    *(x) = FPOF(SIGN(y) | INTOFCONST(INFINITY));

/* Stochastic rounding with equal probabilities. */
#define RS_EQUI_SCALAR(x, y, p, lp)                                            \
  UPDATE_LOCAL_PARAMS(y, p, lp);                                               \
  if (ABS(y) < p->ftzthreshold && *(y) != 0) {         /* Underflow */         \
    randombit = GENBIT(p->BITSEED);                                            \
    *(x) = FPOF(SIGN(y) | INTOFCONST(randombit ? p->ftzthreshold : 0));        \
  } else if (ABS(y) > p->xmax && ABS(y) != INFINITY) { /* Overflow */          \
    randombit = GENBIT(p->BITSEED);                                            \
    *(x) = FPOF(SIGN(y) | INTOFCONST(randombit ? INFINITY : p->xmax));         \
  } else if ((INTOF(y) & lp->trailmask)) { /* Not exactly representable. */    \
    randombit = GENBIT(p->BITSEED);                                            \
    *(x) = FPOF(INTOF(y) & lp->leadmask);                                      \
    if (randombit)                                                             \
      *(x) = FPOF(INTOF(x) + (INTCONST(1) << (DEFPREC-lp->precision)));        \
  } else /* Exactly representable, no rounding necessary. */                   \
    *(x) = *(y);

/* Round-to-odd. */
#define RO_SCALAR(x, y, p, lp)                                                 \
  if (ABS(y) < p->ftzthreshold && *(y) != 0) {         /* Underflow */         \
    *(x) =  FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));                       \
  } else if (ABS(y) > p->xmax && ABS(y) != INFINITY) { /* Overflow */          \
    *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmax));                                \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(y, p, lp);                                             \
    if ((lp->trailmask & INTOF(y)) /* Not exactly representable. */            \
        && (lp->precision > 1)) /* Set last bit to 1 if stored explicitly. */  \
      *(x) = FPOF((INTOF(y) & lp->leadmask) |                                  \
                  (INTCONST(1) << (DEFPREC-lp->precision)));                   \
    else                                                                       \
      *(x) = FPOF(INTOF(y) & lp->leadmask);                                    \
  }

/******************************
 * Macros for vector rounding *
 ******************************/

/*
 * These macros round arrays of type FPTYPE. The input arguments are:
 *   + FPTYPE *X, pointer to memory for storing rounded number;
 *   + FPTYPE *Y, pointer to memory where the number to be rounded is stored;
 *   + (string) PREPROC, pre-processing as *(Y+i) = f(*(A+i), *(B+i), ...);
 *   + (string) POSTPROC, post-processing as *(X+i) = f(*(X+i));
 *   + (string) PARALLEL, either SEQ (for sequential) or PAR (for parallel);
 *   + size_t numlelem, number of elements in the arrays;
 *   + FPPARAMS *p, address where the target format parameters are stored;
 *   + LOCPARAMS *lp, address where the parameters for subnormals are stored.
 *
 * Each element of the array is rounded using one of the scalar rounding macros
 * above. When both a _SAME_EXP and an _OTHER_EXP variant of the scalar rounding
 * function are available, the choice is made by checking whether the maximum
 * exponent of the target format (p->emax) is the same as the maximum exponent
 * of the storage format (DEFEMAX).
 */

#define PARALLEL_STRING(PARALLEL) PARALLEL_STRING_##PARALLEL
#define PARALLEL_STRING_SEQ
#define PARALLEL_STRING_PAR                                                    \
  _Pragma("omp parallel firstprivate(params) private (lp)")

#define FOR_STRING(PARALLEL) FOR_STRING_##PARALLEL
#define FOR_STRING_SEQ
#define FOR_STRING_PAR _Pragma("omp for")

#define VECTOR_ROUNDING_FUNCTION(X, Y, PREPROC, POSTPROC, SCALARFUN, PARALLEL, \
                                 numelem, p, lp)                               \
  PARALLEL_STRING(PARALLEL)                                                    \
  {                                                                            \
    if (p->emax == DEFEMAX) {                                                  \
      FOR_STRING(PARALLEL)                                                     \
      for (size_t i=0; i<numelem; i++) {                                       \
        DEPARENTHESIZE_MAYBE(PREPROC)                                          \
        CONCATENATE(SCALARFUN, _SCALAR_SAME_EXP)(X, Y, p)                      \
        DEPARENTHESIZE_MAYBE(POSTPROC)                                         \
      }                                                                        \
    } else {                                                                   \
      FOR_STRING(PARALLEL)                                                     \
      for (size_t i=0; i<numelem; i++) {                                       \
        DEPARENTHESIZE_MAYBE(PREPROC)                                          \
        CONCATENATE(SCALARFUN, _SCALAR_OTHER_EXP)(X, Y, p, lp)                 \
        DEPARENTHESIZE_MAYBE(POSTPROC)                                         \
      }                                                                        \
    }                                                                          \
  }

/* Round-to-nearest with ties-to-away. */
#define RN_TIES_TO_AWAY(X, Y, PREPROC, POSTPROC, PARALLEL, numelem, p, lp)     \
  VECTOR_ROUNDING_FUNCTION(X, Y, PREPROC, POSTPROC, RN_TIES_TO_AWAY, PARALLEL, \
                           numelem, p, lp)

/* Round-to-nearest with ties-to-zero. */
#define RN_TIES_TO_ZERO(X, Y, PREPROC, POSTPROC, PARALLEL, numelem, p, lp)     \
  VECTOR_ROUNDING_FUNCTION(X, Y, PREPROC, POSTPROC, RN_TIES_TO_ZERO, PARALLEL, \
                           numelem, p, lp)

/* Round-to-nearest with ties-to-even. */
#define RN_TIES_TO_EVEN(X, Y, PREPROC, POSTPROC, PARALLEL, numelem, p, lp) \
  VECTOR_ROUNDING_FUNCTION(X, Y, PREPROC, POSTPROC, RN_TIES_TO_EVEN, PARALLEL, \
                           numelem, p, lp)

/* Round-to-plus-infinity (also known as round-up). */
#define RD_TWD_PINF(X, Y, PREPROC, POSTPROC, PARALLEL, numelem, p, lp)         \
  VECTOR_ROUNDING_FUNCTION(X, Y, PREPROC, POSTPROC, RD_TWD_PINF, PARALLEL,     \
                           numelem, p, lp)

/* Round-to-minus-infinity (also known as round-down). */
#define RD_TWD_NINF(X, Y, PREPROC, POSTPROC, PARALLEL, numelem, p, lp)         \
  VECTOR_ROUNDING_FUNCTION(X, Y, PREPROC, POSTPROC, RD_TWD_NINF, PARALLEL,     \
                           numelem, p, lp)

/* Round-to-zero (also known as truncation). */
#define RD_TWD_ZERO(X, Y, PREPROC, POSTPROC, PARALLEL, numelem, p, lp)         \
  VECTOR_ROUNDING_FUNCTION(X, Y, PREPROC, POSTPROC, RD_TWD_ZERO, PARALLEL,     \
                           numelem, p, lp)

/* Stochastic rounding with proportional probabilities ("mode 1"). */
#define RS_PROP(X, Y, PREPROC, POSTPROC, PARALLEL, INIT_RANDSEED,              \
                numelem, p, lp)                                                \
  PRNG_RAND_INIT                                                               \
  PARALLEL_STRING(PARALLEL)                                                    \
  {                                                                            \
    PRNG_RAND_ADVANCE(PARALLEL)                                                \
    for (size_t i=istart; i<iend; i++) {                                       \
      DEPARENTHESIZE_MAYBE(PREPROC)                                            \
      RS_PROP_SCALAR(X, Y, p, lp)                                              \
      DEPARENTHESIZE_MAYBE(POSTPROC)                                           \
    }                                                                          \
    PRNG_RAND_UPDATE(PARALLEL)                                                 \
  }                                                                            \

/* Stochastic rounding with equal probabilities ("mode 2"). */
#define RS_EQUI(X, Y, PREPROC, POSTPROC, PARALLEL, INIT_BITSEED,               \
                numelem, p, lp)                                                \
  PRNG_BIT_INIT                                                                \
  BITTYPE randombit;                                                           \
  PARALLEL_STRING(PARALLEL)                                                    \
  {                                                                            \
    PRNG_BIT_ADVANCE(PARALLEL)                                                 \
    for (size_t i=istart; i<iend; i++) {                                       \
      DEPARENTHESIZE_MAYBE(PREPROC)                                            \
      RS_EQUI_SCALAR(X, Y, p, lp)                                              \
      DEPARENTHESIZE_MAYBE(POSTPROC)                                           \
    }                                                                          \
    PRNG_BIT_UPDATE(PARALLEL)                                                  \
  }

/* Round-to-odd. */
#define RO(X, Y, PREPROC, POSTPROC, PARALLEL, numelem, p, lp)                  \
  PARALLEL_STRING(PARALLEL)                                                    \
  {                                                                            \
    FOR_STRING(PARALLEL)                                                       \
    for (size_t i=0; i<numelem; i++) {                                         \
      DEPARENTHESIZE_MAYBE(PREPROC)                                            \
      RO_SCALAR(X, Y, p, lp)                                                   \
      DEPARENTHESIZE_MAYBE(POSTPROC)                                           \
    }                                                                          \
  }

/* No rounding. */
#define NO_ROUND(X, Y, PREPROC, POSTPROC, PARALLEL, numelem)                   \
  PARALLEL_STRING(PARALLEL)                                                    \
  {                                                                            \
    FOR_STRING(PARALLEL)                                                       \
    for (size_t i=0; i<numelem; i++) {                                         \
      DEPARENTHESIZE_MAYBE(PREPROC)                                            \
      *(X) = *(Y);                                                             \
      DEPARENTHESIZE_MAYBE(POSTPROC)                                           \
    }                                                                          \
  }

/********************************************************
 * Macros that define rounding or arithmetic functions. *
 ********************************************************/

#define GENERATE_FUN_NAME(FUNNAME)                                             \
  ADDSUFFIXTO(CONCATENATE(CONCATENATE(MAINFUNNAME,_), FUNNAME))
#define GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX)                         \
  CONCATENATE(GENERATE_FUN_NAME(FUNNAME), PARALLEL_SUFFIX)

/*
 * Macros that generates CPFloat functions. They will have name:
 *    <MAINFUNNAME>_<FUNNAME><TYPESUFFIX><PARALLEL_SUFFIX>
 * and their argument list will include:
 *    + (string) FOUTPUT, the output array(s);
 *    + (string) FPINPUT, the input arrays(s);
 *    + size_t numelem, the number of elements in the input and output arrays;
 *    + optstruct *fpopts, address of parameters of target format.
 * The pre-processing and post-processing code should go in the strings:
 *    + PREPROC, featuring Y on the LHS and FPINPUT on the RHS;
 *    + POSTPROC, featuring FPOUTPUT on the LHS and X on the RHS.
 *
 * [Name conventions for input, temporary, and output arrays] The names of the
 * arrays used for rounding are fixed: the rounding routines take the vector Y
 * in input and produce the vector X as output.
 */

#define ADD_BITFLIP_CODE(BITFLIP_CODE) ADD_BITFLIP_CODE_##BITFLIP_CODE
#define ADD_BITFLIP_CODE_NO_BITFLIP(PARALLEL_TYPE, X)
#define ADD_BITFLIP_CODE_INTRODUCE_BITFLIP(PARALLEL_TYPE, X)                   \
  /* Introduce bit flips anywhere in the binary representation. */             \
  if (fpopts->flip != CPFLOAT_NO_SOFTERR) {                                    \
    PRNG_RAND_INIT                                                             \
    size_t flippos;                                                            \
    const size_t totbits = fpopts->flip == CPFLOAT_FP_SOFTERR ?                \
                             NBITS : DEFPREC;                                  \
    FOR_STRING(PARALLEL_TYPE)                                                  \
    for (size_t i = 0; i < numelem; i++) {                                     \
      UPDATE_LOCAL_PARAMS(X+i, &params, &lp);                                  \
      size_t tailbits = DEFPREC - lp.precision;                                \
      if (GENRAND(fpopts->RANDSEED) / (FPTYPE)MAXRAND < fpopts->p) {           \
        flippos = GENRAND(fpopts->RANDSEED) % (totbits - tailbits) + tailbits; \
        *(X) = FPOF(INTOF(X) ^ INTCONST(1) << flippos);                        \
      }                                                                        \
    }                                                                          \
  }

#define GENERATE_INTERFACE_SUBFUNCTIONS(FUNNAME,                               \
                                        FOUTPUT,                               \
                                        FINPUT,                                \
                                        X, Y,                                  \
                                        INIT_STRING,                           \
                                        BITFLIP_CODE,                          \
                                        PREPROC,                               \
                                        POSTPROC,                              \
                                        PARALLEL_TYPE,                         \
                                        PARALLEL_SUFFIX,                       \
                                        INIT_RAND_STRING,                      \
                                        INIT_BIT_STRING)                       \
  static inline int                                                            \
  GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX)                               \
  (DEPARENTHESIZE_MAYBE(FOUTPUT),                                              \
   DEPARENTHESIZE_MAYBE(FINPUT),                                               \
   const size_t numelem, optstruct *fpopts) {                                  \
    int retval = 0;                                                            \
    FPPARAMS params = COMPUTE_GLOBAL_PARAMS(fpopts, &retval);                  \
    LOCPARAMS lp;                                                              \
    DEPARENTHESIZE_MAYBE(INIT_STRING)                                          \
    switch (params.round) {                                                    \
    case CPFLOAT_RND_NA: /* round-to-nearest with ties-to-away  */             \
      RN_TIES_TO_AWAY(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE,                  \
                      numelem, (&params), (&lp));                              \
      break;                                                                   \
    case CPFLOAT_RND_NZ: /* round-to-nearest with ties-to-zero */              \
      RN_TIES_TO_ZERO(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE,                  \
                      numelem, (&params), (&lp));                              \
      break;                                                                   \
    case CPFLOAT_RND_NE: /* round-to-nearest with ties-to-even */              \
      RN_TIES_TO_EVEN(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE,                  \
                      numelem, (&params), (&lp));                              \
      break;                                                                   \
    case CPFLOAT_RND_TP: /* round-toward-positive */                           \
      RD_TWD_PINF(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE,                      \
                      numelem, (&params), (&lp));                              \
      break;                                                                   \
    case CPFLOAT_RND_TN: /* round-toward-negative */                           \
      RD_TWD_NINF(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE,                      \
                      numelem, (&params), (&lp));                              \
      break;                                                                   \
    case CPFLOAT_RND_TZ: /* round-toward-zero */                               \
      RD_TWD_ZERO(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE,                      \
                      numelem, (&params), (&lp));                              \
      break;                                                                   \
    case CPFLOAT_RND_SP: /* stochastic rounding (proportional) */              \
      RS_PROP(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE, INIT_RAND_STRING,        \
              numelem, (&params), (&lp));                                      \
      break;                                                                   \
    case CPFLOAT_RND_SE: /* stochastic rounding (equal) */                     \
      RS_EQUI(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE, INIT_BIT_STRING,         \
              numelem, (&params), (&lp));                                      \
      break;                                                                   \
    case CPFLOAT_RND_OD: /* round-to-odd */                                    \
      RO(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE,                               \
         numelem, (&params), (&lp));                                           \
      break;                                                                   \
    default: /* No rounding if unknown mode specified. */                      \
      NO_ROUND(X, Y, PREPROC, POSTPROC, PARALLEL_TYPE, numelem);               \
      break;                                                                   \
    }                                                                          \
    ADD_BITFLIP_CODE(BITFLIP_CODE)(PARALLEL_TYPE, X)                           \
    return retval;                                                             \
  }

#define PARALLEL_SUFFIX_SEQ _sequential
#define INIT_RAND_STRING_SEQ(p) INIT_RANDSEED_SEQ(p);
#define INIT_BIT_STRING_SEQ(p) INIT_BITSEED_SEQ(p);

#define GENERATE_INTERFACE_SUBFUNCTIONS_SEQ(FUNNAME,                           \
                                            FOUTPUT,                           \
                                            FINPUT,                            \
                                            X, Y,                              \
                                            INIT_STRING, FINAL_STRING,         \
                                            PREPROC, POSTPROC)                 \
  GENERATE_INTERFACE_SUBFUNCTIONS(FUNNAME,                                     \
                                  FOUTPUT,                                     \
                                  FINPUT,                                      \
                                  X, Y,                                        \
                                  INIT_STRING,                                 \
                                  FINAL_STRING,                                \
                                  PREPROC,                                     \
                                  POSTPROC,                                    \
                                  SEQ,                                         \
                                  PARALLEL_SUFFIX_SEQ,                         \
                                  INIT_RAND_STRING_SEQ,                        \
                                  INIT_BIT_STRING_SEQ)

#ifdef _OPENMP
#define PARALLEL_SUFFIX_PAR _parallel
#define INIT_RAND_STRING_PAR INIT_RANDSEED_PAR
#define INIT_BIT_STRING_PAR INIT_BITSEED_PAR

#define GENERATE_INTERFACE_SUBFUNCTIONS_PAR(FUNNAME,                           \
                                            FOUTPUT,                           \
                                            FINPUT,                            \
                                            X, Y,                              \
                                            INIT_STRING, FINAL_STRING,         \
                                            PREPROC, POSTPROC)                 \
  GENERATE_INTERFACE_SUBFUNCTIONS(FUNNAME,                                     \
                                  FOUTPUT,                                     \
                                  FINPUT,                                      \
                                  X, Y,                                        \
                                  INIT_STRING,                                 \
                                  FINAL_STRING,                                \
                                  PREPROC,                                     \
                                  POSTPROC,                                    \
                                  PAR,                                         \
                                  PARALLEL_SUFFIX_PAR,                         \
                                  INIT_RAND_STRING_PAR,                        \
                                  INIT_BIT_STRING_PAR)

/* Partial results before rounding are stored in the vector X. */
#define GENERATE_INTERFACE(FUNNAME, FOUTPUT, FINPUT, X, Y,                     \
                           FCALLOUT, FCALLIN,                                  \
                           INIT_STRING, FINAL_STRING, PREPROC, POSTPROC)       \
  GENERATE_INTERFACE_SUBFUNCTIONS_SEQ(FUNNAME, FOUTPUT, FINPUT, X, Y,          \
                                      INIT_STRING, FINAL_STRING,               \
                                      PREPROC, POSTPROC)                       \
  GENERATE_INTERFACE_SUBFUNCTIONS_PAR(FUNNAME, FOUTPUT, FINPUT, X, Y,          \
                                      INIT_STRING, FINAL_STRING,               \
                                      PREPROC, POSTPROC)                       \
  static inline int                                                            \
  GENERATE_FUN_NAME(FUNNAME)(DEPARENTHESIZE_MAYBE(FOUTPUT),                    \
                             DEPARENTHESIZE_MAYBE(FINPUT),                     \
                             const size_t numelem,                             \
                             optstruct *fpopts) {                              \
    if (numelem < CONCATENATE(OPENMP_THRESHOLD_, FPTYPE))                      \
      return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_SEQ)                \
        (DEPARENTHESIZE_MAYBE(FCALLOUT),                                       \
         DEPARENTHESIZE_MAYBE(FCALLIN),                                        \
         numelem, fpopts);                                                     \
    else                                                                       \
      return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_PAR)                \
        (DEPARENTHESIZE_MAYBE(FCALLOUT),                                       \
         DEPARENTHESIZE_MAYBE(FCALLIN),                                        \
         numelem, fpopts);                                                     \
  }
#else /* #ifdef _OPENMP */
#define GENERATE_INTERFACE(FUNNAME, FOUTPUT, FINPUT, X, Y,                     \
                           FCALLOUT, FCALLIN,                                  \
                           INIT_STRING, FINAL_STRING,                          \
                           PREPROC, POSTPROC)                                  \
  GENERATE_INTERFACE_SUBFUNCTIONS_SEQ(FUNNAME, FOUTPUT, FINPUT, X, Y,          \
                                      INIT_STRING, FINAL_STRING,               \
                                      PREPROC, POSTPROC)                       \
  static inline int                                                            \
  GENERATE_FUN_NAME(FUNNAME)(DEPARENTHESIZE_MAYBE(FOUTPUT),                    \
                             DEPARENTHESIZE_MAYBE(FINPUT),                     \
                             const size_t numelem,                             \
                             optstruct *fpopts) {                              \
    return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_SEQ)                  \
      (DEPARENTHESIZE_MAYBE(FCALLOUT),                                         \
       DEPARENTHESIZE_MAYBE(FCALLIN),                                          \
       numelem, fpopts);                                                       \
  }
#endif /* #ifdef _OPENMP */

/* Function generators.
 *
 * If necessary, the generated functions use the output matrix X to store the
 * intermediate results of the computation, then round in place.
*/
#define GENERATE_UNIVARIATE(FUNNAME, OPSTRING)                                 \
  GENERATE_INTERFACE(FUNNAME, FPTYPE *X,                                       \
                     const FPTYPE *A,                                          \
                     X+i, X+i,                                                 \
                     X, A,                                                     \
                     NOOP, INTRODUCE_BITFLIP,                                  \
                     *(X+i) = OPSTRING(*(A+i));, NOOP)

#define GENERATE_BIVARIATE_PREFIX(FUNNAME, OPSTRING)                           \
  GENERATE_INTERFACE(FUNNAME, FPTYPE *X,                                       \
                     (const FPTYPE *A, const FPTYPE *B),                       \
                     X+i, X+i,                                                 \
                     X, (A, B),                                                \
                     NOOP, INTRODUCE_BITFLIP,                                  \
                     *(X+i) = OPSTRING(*(A+i), *(B+i));, NOOP)

#define GENERATE_BIVARIATE_INFIX(FUNNAME, OPSTRING)                            \
  GENERATE_INTERFACE(FUNNAME, FPTYPE *X,                                       \
                     (const FPTYPE *A, const FPTYPE *B),                       \
                     X+i, X+i,                                                 \
                     X, (A, B),                                                \
                     NOOP, INTRODUCE_BITFLIP,                                  \
                     *(X+i) = *(A+i) OPSTRING *(B+i);, NOOP)

#define GENERATE_TRIVARIATE(FUNNAME, OPSTRING)                                 \
  GENERATE_INTERFACE(FUNNAME, FPTYPE *X,                                       \
                     (const FPTYPE *A, const FPTYPE *B, const FPTYPE *C),      \
                     X+i, X+i,                                                 \
                     X, (A, B, C),                                             \
                     NOOP, INTRODUCE_BITFLIP,                                  \
                     *(X+i) = OPSTRING(*(A+i), *(B+i), *(C+i));, NOOP)

#define GENERATE_UNIVARIATE_POSTPROCESSING(FUNNAME, OPSTRING)                  \
  GENERATE_INTERFACE(FUNNAME, FPTYPE *X,                                       \
                     const FPTYPE *A,                                          \
                     X+i, A+i,                                                 \
                     X, A,                                                     \
                     NOOP, INTRODUCE_BITFLIP,                                  \
                     NOOP, *(X+i) = OPSTRING(*(X+i));)

/* The following macros use the default names in the math.h library, that is,
   they add the type prefix (f for float, nothing for double) to the name of the
   prefix operator. */
#define GENERATE_UNIVARIATE_MATH_H(FUNNAME, OPSTRING)                          \
  GENERATE_UNIVARIATE(FUNNAME, ADDSUFFIXTO(OPSTRING))
#define GENERATE_BIVARIATE_PREFIX_MATH_H(FUNNAME, OPSTRING)                    \
  GENERATE_BIVARIATE_PREFIX(FUNNAME, ADDSUFFIXTO(OPSTRING))
#define GENERATE_TRIVARIATE_MATH_H(FUNNAME, OPSTRING)                          \
  GENERATE_TRIVARIATE(FUNNAME, ADDSUFFIXTO(OPSTRING))
#define GENERATE_UNIVARIATE_POSTPROCESSING_MATH_H(FUNNAME, OPSTRING)           \
  GENERATE_UNIVARIATE_POSTPROCESSING(FUNNAME, ADDSUFFIXTO(OPSTRING))

/* The following are interfaces for the nextafter and nexttoward functions. */
#define GENERATE_NEXTAFTER_SUBFUNCTIONS(FUNNAME,                               \
                                        TYPETO,                                \
                                        PARALLEL_TYPE,                         \
                                        PARALLEL_SUFFIX)                       \
  static inline int                                                            \
  GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX)                               \
  (FPTYPE *X, const FPTYPE *A, const TYPETO *B,                                \
   const size_t numelem, optstruct *fpopts) {                                  \
    int retval = 0;                                                            \
    FPPARAMS params = COMPUTE_GLOBAL_PARAMS(fpopts, &retval);                  \
    LOCPARAMS lp;                                                              \
    PARALLEL_STRING(PARALLEL_TYPE)                                             \
      {                                                                        \
        FOR_STRING(PARALLEL_TYPE)                                              \
        for (size_t i=0; i<numelem; i++) {                                     \
          int rounddir = *(B+i) > *(A+i);                                      \
          *(X+i) = ADDSUFFIXTO(FUNNAME)(*(A+i), *(B+i));                       \
          if (rounddir) {                                                      \
            RD_TWD_PINF_SCALAR_OTHER_EXP(X+i, X+i, (&params), (&lp))           \
          } else {                                                             \
            RD_TWD_NINF_SCALAR_OTHER_EXP(X+i, X+i, (&params), (&lp))           \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    return retval;                                                             \
  }

#ifdef _OPENMP
#define GENERATE_NEXTAFTER_INTERFACE(FUNNAME, TYPETO)                          \
  GENERATE_NEXTAFTER_SUBFUNCTIONS(FUNNAME, TYPETO, SEQ, PARALLEL_SUFFIX_SEQ)   \
  GENERATE_NEXTAFTER_SUBFUNCTIONS(FUNNAME, TYPETO, PAR, PARALLEL_SUFFIX_PAR)   \
  static inline int                                                            \
       GENERATE_FUN_NAME(FUNNAME)(FPTYPE *X, const FPTYPE *A, const TYPETO *B, \
                             const size_t numelem, optstruct *fpopts) {        \
    if (numelem < CONCATENATE(OPENMP_THRESHOLD_, FPTYPE))                      \
      return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_SEQ)                \
        (X, A, B, numelem, fpopts);                                            \
    else                                                                       \
      return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_PAR)                \
        (X, A, B, numelem, fpopts);                                            \
  }
#else /* #ifdef _OPENMP */
#define GENERATE_NEXTAFTER_INTERFACE(FUNNAME, TYPETO)                          \
  GENERATE_NEXTAFTER_SUBFUNCTIONS(FUNNAME, TYPETO, SEQ, PARALLEL_SUFFIX_SEQ)   \
  static inline int                                                            \
  GENERATE_FUN_NAME(FUNNAME)(FPTYPE *X, const FPTYPE *A, const TYPETO *B,      \
                             const size_t numelem, optstruct *fpopts) {        \
      return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_SEQ)                \
        (X, A, B, numelem, fpopts);                                            \
  }
#endif /* #ifdef _OPENMP */

/**************************************
 * Actual instantiation of functions. *
 **************************************/

/* Only rounding. */
GENERATE_INTERFACE(fpround,
                   FPTYPE *X,
                   const FPTYPE *A,
                   X+i, A+i,
                   X, A,
                   NOOP, INTRODUCE_BITFLIP, NOOP, NOOP)
static inline
int ADDSUFFIXTO(cpfloat)(FPTYPE *X,
                         const FPTYPE *A,
                         const size_t numelem,
                         optstruct *fpopts) {
  return GENERATE_FUN_NAME(fpround)(X, A, numelem, fpopts);
}

/*
 * The following two functions are used internally to autotune the size of the
 * smallest vector on which multiple OpenMP threads can be used, and to perform
 * some of the experiments in the experiments/ folder.
 * They are not documented and should not be used directly.
 */
#ifdef _OPENMP
static inline
int CONCATENATE(ADDSUFFIXTO(cpfloat),_sequential)(FPTYPE *X,
                                                  const FPTYPE *A,
                                                  const size_t numelem,
                                                  optstruct *fpopts) {
  return GENERATE_SUBFUN_NAME(fpround, PARALLEL_SUFFIX_SEQ)
    (X, A, numelem, fpopts);
}
static inline
int CONCATENATE(ADDSUFFIXTO(cpfloat),_parallel)(FPTYPE *X,
                                                const FPTYPE *A,
                                                const size_t numelem,
                                                optstruct *fpopts) {
  return GENERATE_SUBFUN_NAME(fpround, PARALLEL_SUFFIX_PAR)
    (X, A, numelem, fpopts);
}
#endif /* #ifdef _OPENMP */

/* Elementary arithmetic operations. */
GENERATE_BIVARIATE_INFIX(add, +)
GENERATE_BIVARIATE_INFIX(sub, -)
GENERATE_BIVARIATE_INFIX(mul, *)
GENERATE_BIVARIATE_INFIX(div, /)

/* Trigonometric functions. */
GENERATE_UNIVARIATE_MATH_H(cos, cos)
GENERATE_UNIVARIATE_MATH_H(sin, sin)
GENERATE_UNIVARIATE_MATH_H(tan, tan)

GENERATE_UNIVARIATE_MATH_H(acos, acos)
GENERATE_UNIVARIATE_MATH_H(asin, asin)
GENERATE_UNIVARIATE_MATH_H(atan, atan)
GENERATE_BIVARIATE_PREFIX(atan2, atan2)

/* Hyperbolic functions. */
GENERATE_UNIVARIATE_MATH_H(cosh, cosh)
GENERATE_UNIVARIATE_MATH_H(sinh, sinh)
GENERATE_UNIVARIATE_MATH_H(tanh, tanh)

GENERATE_UNIVARIATE_MATH_H(acosh, acosh)
GENERATE_UNIVARIATE_MATH_H(asinh, asinh)
GENERATE_UNIVARIATE_MATH_H(atanh, atanh)

/* Exponentiation and logarithmic functions. */
GENERATE_UNIVARIATE_MATH_H(exp, exp)
GENERATE_INTERFACE(frexp,
                   (FPTYPE *X, int *exp), const FPTYPE *A,
                   X+i, A+i,
                   (X, exp), A,
                   NOOP, INTRODUCE_BITFLIP,
                   NOOP, *(X+i) = ADDSUFFIXTO(frexp)(*(X+i), exp+i);)
GENERATE_INTERFACE(ldexp,
                   FPTYPE *X, (const FPTYPE *A, const int *exp),
                   X+i, X+i,
                   X, (A, exp),
                   NOOP, INTRODUCE_BITFLIP,
                   *(X+i) = ADDSUFFIXTO(ldexp)(*(A+i), *(exp+i));, NOOP)
GENERATE_UNIVARIATE_MATH_H(log, log)
GENERATE_UNIVARIATE_MATH_H(log10, log10)
GENERATE_INTERFACE(modf,
                   (FPTYPE *X, FPTYPE *intpart), const FPTYPE *A,
                   X+i, A+i,
                   (X, intpart), A,
                   NOOP, INTRODUCE_BITFLIP,
                   NOOP, *(X+i) = ADDSUFFIXTO(modf)(*(X+i), intpart+i);)
GENERATE_UNIVARIATE_MATH_H(exp2, exp2)
GENERATE_UNIVARIATE_MATH_H(expm1, expm1)
GENERATE_INTERFACE(ilogb,
                   int *exp, const FPTYPE *A,
                   &loctmp, A+i,
                   exp, A,
                   NOOP, NO_BITFLIP,
                   FPTYPE loctmp;, *(exp+i) =  ADDSUFFIXTO(ilogb)(loctmp);)

GENERATE_UNIVARIATE_MATH_H(log1p, log1p)
GENERATE_UNIVARIATE_MATH_H(log2, log2)
GENERATE_INTERFACE(logb,
                   FPTYPE *X, const FPTYPE *A,
                   X+i, A+i,
                   X, A,
                   NOOP, INTRODUCE_BITFLIP,
                   NOOP, *(X+i) = ADDSUFFIXTO(logb)(*(X+i));)
GENERATE_INTERFACE(scalbn,
                   FPTYPE *X, (const FPTYPE *A, const int *exp),
                   X+i, A+i,
                   X, (A, exp),
                   NOOP, INTRODUCE_BITFLIP,
                   NOOP, *(X+i) = ADDSUFFIXTO(scalbn)(*(X+i), *(exp+i));)
GENERATE_INTERFACE(scalbln,
                   FPTYPE *X, (const FPTYPE *A, const long int *exp),
                   X+i, A+i,
                   X, (A, exp),
                   NOOP, INTRODUCE_BITFLIP,
                   NOOP, *(X+i) = ADDSUFFIXTO(scalbln)(*(X+i), *(exp+i));)

/* Power functions. */
GENERATE_BIVARIATE_PREFIX_MATH_H(pow, pow)
GENERATE_UNIVARIATE_MATH_H(sqrt, sqrt)
GENERATE_UNIVARIATE_MATH_H(cbrt, cbrt)
GENERATE_BIVARIATE_PREFIX_MATH_H(hypot, hypot)

/* Error and gamma functions. */
GENERATE_UNIVARIATE_MATH_H(erf, erf)
GENERATE_UNIVARIATE_MATH_H(erfc, erfc)
GENERATE_UNIVARIATE_MATH_H(tgamma, tgamma)
GENERATE_UNIVARIATE_MATH_H(lgamma, lgamma)

/* Rounding and remainder functions */
#define FIND_ROUNDING_PRECISION(extra_instructions)                            \
  (extra_instructions                                                          \
   int temp = ADDSUFFIXTO(ilogb)(*(A+i));                                      \
   if (temp == FP_ILOGB0 || temp == FP_ILOGBNAN || temp == INT_MAX) {          \
     params.precision = DEFPREC-1;                                             \
     params.emax = DEFEMAX;                                                    \
     params.emax = DEFEMIN;                                                    \
   } else {                                                                    \
     params.precision = temp + 1;                                              \
     params.ftzthreshold = 1.0;                                                \
     params.xmin = 1.0;                                                        \
   }                                                                           \
   UPDATE_GLOBAL_BITMASKS((&params));)

#define RINT_INIT_STRING                                                       \
  (if (fpopts->round == CPFLOAT_RND_NE) {                                      \
    params.subnormal = CPFLOAT_SUBN_USE;                                       \
    params.emax = DEFEMAX - 1;                                                 \
  } else {                                                                     \
    params.subnormal = CPFLOAT_SUBN_RND;                                       \
    params.emax = DEFEMAX;                                                     \
  }                                                                            \
  params.xmax = ldexp(1., fpopts->emax) * (2-ldexp(1., 1-fpopts->precision));)

GENERATE_INTERFACE(ceil,
                   FPTYPE *X, const FPTYPE *A,
                   X+i, A+i,
                   X, A,
                   params.round = CPFLOAT_RND_TP;, INTRODUCE_BITFLIP,
                   FIND_ROUNDING_PRECISION(),
                   NOOP)
GENERATE_INTERFACE(floor,
                   FPTYPE *X, const FPTYPE *A,
                   X+i, A+i,
                   X, A,
                   params.round = CPFLOAT_RND_TN;, INTRODUCE_BITFLIP,
                   FIND_ROUNDING_PRECISION(),
                   NOOP)
GENERATE_BIVARIATE_PREFIX_MATH_H(fmod, fmod)
GENERATE_INTERFACE(trunc,
                   FPTYPE *X, const FPTYPE *A,
                   X+i, A+i,
                   X, A,
                   params.round = CPFLOAT_RND_TZ;, INTRODUCE_BITFLIP,
                   FIND_ROUNDING_PRECISION(),
                   NOOP)

GENERATE_INTERFACE(round,
                   FPTYPE *X, const FPTYPE *A,
                   X+i, A+i,
                   X, A,
                   params.round = CPFLOAT_RND_NA;, INTRODUCE_BITFLIP,
                   FIND_ROUNDING_PRECISION(),
                   NOOP)
GENERATE_INTERFACE(lround,
                   long *r, const FPTYPE *A,
                   &loctmp, A+i,
                   r, A,
                   params.round = CPFLOAT_RND_NA;, NO_BITFLIP,
                   FIND_ROUNDING_PRECISION(FPTYPE loctmp;),
                   *(r+i) = (long)loctmp;)
GENERATE_INTERFACE(llround,
                   long long *r, const FPTYPE *A,
                   &loctmp, A+i,
                   r, A,
                   params.round = CPFLOAT_RND_NA;, NO_BITFLIP,
                   FIND_ROUNDING_PRECISION(FPTYPE loctmp;),
                   *(r+i) = (long long)loctmp;)

GENERATE_INTERFACE(rint,
                   (FPTYPE *X, int *exception), const FPTYPE *A,
                   X+i, A+i,
                   (X, exception), A,
                   RINT_INIT_STRING, INTRODUCE_BITFLIP,
                   FIND_ROUNDING_PRECISION(),
                   *(exception+i) = *(X+i) == *(A+i) ? 0 : FE_INEXACT;)
GENERATE_INTERFACE(lrint,
                   (long *r, int *exception),
                   const FPTYPE *A,
                   &loctmp, A+i,
                   (r, exception), A,
                   RINT_INIT_STRING, NO_BITFLIP,
                   FIND_ROUNDING_PRECISION(FPTYPE loctmp;),
                   (*(r+i) = (long)loctmp;
                    *(exception+i) = *(r+i) == *(A+i) ? 0 : FE_INEXACT;))
GENERATE_INTERFACE(llrint,
                   (long long *r, int *exception),
                   const FPTYPE *A,
                   &loctmp, A+i,
                   (r, exception), A,
                   RINT_INIT_STRING, NO_BITFLIP,
                   FIND_ROUNDING_PRECISION(FPTYPE loctmp;),
                   (*(r+i) = (long long)loctmp;
                    *(exception+i) = *(r+i) == *(A+i) ? 0 : FE_INEXACT;))
GENERATE_INTERFACE(nearbyint,
                   FPTYPE *X, const FPTYPE *A,
                   X+i, A+i,
                   X, A,
                   RINT_INIT_STRING, INTRODUCE_BITFLIP,
                   FIND_ROUNDING_PRECISION(),
                   NOOP)
GENERATE_BIVARIATE_PREFIX_MATH_H(remainder, remainder)

GENERATE_INTERFACE(remquo,
                   (FPTYPE *X, int *quot),
                   (const FPTYPE *A, const FPTYPE *B),
                   X+i, X+i,
                   (X, quot), (A, B),
                   NOOP, INTRODUCE_BITFLIP,
                   *(X+i) = ADDSUFFIXTO(remquo)(*(A+i), *(B+i), quot+i);,
                   NOOP)

/* Floating-point manipulation functions. */
GENERATE_BIVARIATE_PREFIX_MATH_H(copysign, copysign)
/* nan NOT IMPLEMENTED as not relevant. */
GENERATE_NEXTAFTER_INTERFACE(nextafter, FPTYPE)
GENERATE_NEXTAFTER_INTERFACE(nexttoward, long double)

/* Minimum, maximum, difference functions. */
GENERATE_BIVARIATE_PREFIX_MATH_H(fdim, fdim)
GENERATE_BIVARIATE_PREFIX_MATH_H(fmax, fmax)
GENERATE_BIVARIATE_PREFIX_MATH_H(fmin, fmin)

/* Classification. */
GENERATE_INTERFACE(fpclassify,
                   int *r, const FPTYPE *A,
                   &loctmp, A+i,
                   r, A,
                   NOOP, NO_BITFLIP, FPTYPE loctmp;,
                   if (loctmp == 0)
                     *(r+i) = FP_ZERO;
                   else if (ABS(&loctmp) < params.xmin)
                     *(r+i) = FP_SUBNORMAL;
                   else if (ABS(&loctmp) < INFINITY)
                     *(r+i) = FP_NORMAL;
                   else if (isinf(loctmp))
                     *(r+i) = FP_INFINITE;
                   else
                     *(r+i) = FP_NAN;)
GENERATE_INTERFACE(isfinite,
                   int *r, const FPTYPE *A,
                   &loctmp, A+i,
                   r, A,
                   NOOP, NO_BITFLIP,
                   FPTYPE loctmp;, *(r+i) = isfinite(loctmp);)
GENERATE_INTERFACE(isinf,
                   int *r, const FPTYPE *A,
                   &loctmp, A+i,
                   r, A,
                   NOOP,  NO_BITFLIP,
                   FPTYPE loctmp;, *(r+i) = isinf(loctmp);)
GENERATE_INTERFACE(isnan,
                   int *r, const FPTYPE *A,
                   &loctmp, A+i,
                   r, A,
                   NOOP,  NO_BITFLIP,
                   FPTYPE loctmp;, *(r+i) = isnan(loctmp);)
GENERATE_INTERFACE(isnormal,
                   int *r, const FPTYPE *A,
                   &loctmp, A+i,
                   r, A,
                   NOOP,  NO_BITFLIP, FPTYPE loctmp;,
                   *(r+i) = ABS(&loctmp) >= params.xmin && !isinf(loctmp);)
/* signbit NOT IMPLEMENTED as rounding doesn't interfere with it. */

/* Comparison.
 * Not implemented as they don't seem to be meaningful for vectors.
 *isgreater NOT IMPLEMENTED
 *isgreaterequal NOT IMPLEMENTED
 *isless NOT IMPLEMENTED
 *islessequal NOT IMPLEMENTED
 *islessgreater NOT IMPLEMENTED
 *isunordered NOT IMPLEMENTED
 */

/* Other functions. */
GENERATE_UNIVARIATE_MATH_H(fabs, fabs)
GENERATE_TRIVARIATE_MATH_H(fma, fma)





/* Undefine local macros. */
#undef ADDSUFFIXTO
#undef FPUNION
#undef INTCONST
#undef INTOF
#undef INTOFCONST
#undef FPOF

#undef RANDSEED
#undef RANDSEEDTYPE

#undef FPPARAMS
#undef LOCPARAMS
#undef SIGN
#undef ABS

#ifndef PCG_VARIANTS_H_INCLUDED
#ifdef _OPENMP
#undef PRNG_ADVANCE_RAND
#undef ADVANCERAND
#endif /* #ifdef _OPENMP */
#undef GENRAND
#endif  /* #ifdef PCG_VARIANTS_H_INCLUDED */

#undef PRNG_RAND_INIT

#undef PRNG_RAND_ADVANCE
#undef PRNG_RAND_ADVANCE_SEQ
#undef PRNG_RAND_ADVANCE_PAR

#undef PRNG_RAND_UPDATE
#undef PRNG_RAND_UPDATE_SEQ
#undef PRNG_RAND_UPDATE_PAR


#undef VALIDATE_INPUT
#undef COMPUTE_GLOBAL_PARAMS
#undef UPDATE_LOCAL_PARAMS

#undef RN_TIES_TO_AWAY_SCALAR_SAME_EXP
#undef RN_TIES_TO_AWAY_SCALAR_OTHER_EXP
#undef RN_TIES_TO_ZERO_SCALAR_SAME_EXP
#undef RN_TIES_TO_ZERO_SCALAR_OTHER_EXP
#undef RN_TIES_TO_EVEN_SCALAR_SAME_EXP
#undef RN_TIES_TO_EVEN_SCALAR_OTHER_EXP
#undef RD_TWD_PINF_SCALAR_SAME_EXP
#undef RD_TWD_PINF_SCALAR_OTHER_EXP
#undef RD_TWD_NINF_SCALAR_SAME_EXP
#undef RD_TWD_NINF_SCALAR_OTHER_EXP
#undef RD_TWD_ZERO_SCALAR_SAME_EXP
#undef RD_TWD_ZERO_SCALAR_OTHER_EXP
#undef RS_PROP_SCALAR
#undef RS_EQUI_SCALAR
#undef RO_SCALAR

#undef FOR_STRING
#undef FOR_STRING_SEQ
#undef FOR_STRING_PAR

#undef RN_TIES_TO_AWAY
#undef RN_TIES_TO_ZERO
#undef RN_TIES_TO_EVEN
#undef RD_TWD_PINF
#undef RD_TWD_NINF
#undef RD_TWD_ZERO
#undef RS_PROP
#undef RS_EQUI
#undef RO

#undef GENERATE_FUN_NAME
#undef GENERATE_SUBFUN_NAME
#undef PARALLEL_STRING
#undef PARALLEL_STRING_SEQ
#undef PARALLEL_STRING_PAR

#undef GENERATE_INTERFACE_SUBFUNCTIONS
#undef PARALLEL_SUFFIX_SEQ
#undef INIT_RAND_STRING_SEQ
#undef INIT_BIT_STRING_SEQ
#undef GENERATE_INTERFACE_SUBFUNCTIONS_SEQ

#undef PARALLEL_SUFFIX_PAR
#undef INIT_RAND_STRING_PAR
#undef INIT_BIT_STRING_PAR
#undef GENERATE_INTERFACE_SUBFUNCTIONS_PAR

#undef GENERATE_INTERFACE
#undef GENERATE_UNIVARIATE
#undef GENERATE_BIVARIATE_PREFIX
#undef GENERATE_BIVARIATE_INFIX
#undef GENERATE_TRIVARIATE
#undef GENERATE_UNIVARIATE_POSTPROCESSING

#undef GENERATE_UNIVARIATE_MATH_H
#undef GENERATE_BIVARIATE_PREFIX_MATH_H
#undef GENERATE_TRIVARIATE_MATH_H
#undef GENERATE_UNIVARIATE_POSTPROCESSING_MATH_H

/* Undefine macros from external file. */
#undef FUNSUFFIX
#undef FPTYPE
#undef INTTYPE
#undef INTSUFFIX

#undef DEFPREC
#undef DEFEMAX
#undef DEFEMIN
#undef NLEADBITS
#undef NBITS
#undef FULLMASK
#undef ABSMASK
#undef SIGNMASK
#undef EXPMASK
#undef FRACMASK

#undef MAXRAND
#undef INITRAND
#ifdef PCG_VARIANTS_H_INCLUDED
#undef ADVANCERAND
#undef GENRAND
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#undef GEN_SINGLE_RAND
#endif /* #ifdef PCG_VARIANTS_H_INCLUDED */

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
