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

/****************************
 ****************************
 * TYPE-INDEPENDENT MACROS. *
 ****************************
 ****************************/

/* Define a no-operation empty macro. */
#define NOOP

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
#define INITBIT_SEQ(seed)                       \
  pcg32_srandom_r(seed,                         \
                  time(NULL),                   \
                  (intptr_t)seed)
#ifdef _OPENMP
#define INITBIT_PAR(seed)                                       \
  pcg32_srandom_r(seed,                                         \
                  omp_get_thread_num() * 13254 + time(NULL),    \
                  (intptr_t)seed)
#endif /* #ifdef _OPENMP */
#define GENBIT(seed) (pcg32_random_r(seed) & (1U << 31))
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#define INITBIT_SEQ(seed) *seed = time(NULL)
#ifdef _OPENMP
#define INITBIT_PAR(seed) *seed = omp_get_thread_num() * 13254 + time(NULL)
#endif /* #ifdef _OPENMP */
#define GENBIT(seed) (rand_r(seed) & (1U << 30))
#endif /* #ifdef PCG_VARIANTS_H_INCLUDED */

/* Prefix of all library functions. */
#define MAINFUNNAME cpfloat





/**************************
 **************************
 * TYPE-DEPENDENT MACROS. *
 **************************
 **************************/

#define ADDSUFFIXTO(x) CONCATENATE(x, FUNSUFFIX)

/* Struct to reinterpret floating-point values as unsigned integers. */
#define FPUNION ADDSUFFIXTO(fpint)
typedef union {
  INTTYPE intval;
  FPTYPE fpval;
} FPUNION;

/* Types for internal state of pseudo-random number generator. */
#define RANDSEED ADDSUFFIXTO(randseed)
#define RANDSEEDTYPE CONCATENATE(ADDSUFFIXTO(cpfloat_randseed),_t)

/* Relevant parameters of floating-point number system. */
#define FPPARAMS ADDSUFFIXTO(fpparams)
typedef struct {
  const cpfloat_precision_t precision;
  const cpfloat_exponent_t emax;
  const cpfloat_exponent_t emin;
  const cpfloat_subnormal_t subnormal;
  const FPTYPE ftzthreshold;
  const FPTYPE xmin;
  const FPTYPE xmax;
  const FPTYPE xbnd;
  const INTTYPE leadmask;
  const INTTYPE trailmask;
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

#define INTCONST(x) CONCATENATE(x, INTSUFFIX)

#define INTOF(x)(((FPUNION *)(x))->intval)
#define INTOFCONST(x)(((FPUNION){.fpval = (FPTYPE)x}).intval)
#define FPOF(x)(((FPUNION){.intval = (INTTYPE)(x)}).fpval)

#define SIGN(x)(SIGNMASK & INTOF(x))
#define ABS(x)(FPOF((INTTYPE)(ABSMASK & INTOF(x))))

/* Functions to initialize the pseudo-random number generator. */
#ifdef _OPENMP
#define INIT_RANDSEED_PAR CONCATENATE(ADDSUFFIXTO(init_randseed), _par)
static inline void INIT_RANDSEED_PAR (FPPARAMS *params) {
  if (params->RANDSEED == NULL) {
    params->RANDSEED = malloc(sizeof(*params->RANDSEED));
    INITRAND_PAR(params->RANDSEED);
  }
}
#endif /* #ifdef _OPENMP */

/* Functions to initialize the pseudo-random bit generator. */
#ifdef _OPENMP
#define INIT_BITSEED_PAR CONCATENATE(ADDSUFFIXTO(init_bitseed), _par)
static inline void INIT_BITSEED_PAR (FPPARAMS *params) {
  if (params->BITSEED == NULL) {
    params->BITSEED = malloc(sizeof(*params->BITSEED));
    INITBIT_PAR(params->BITSEED);
  }
}
#endif /* #ifdef _OPENMP */

#define INIT_BITSEED_SEQ CONCATENATE(ADDSUFFIXTO(init_bitseed), _seq)
static inline void INIT_BITSEED_SEQ (FPPARAMS *params) {
  if (params->BITSEED == NULL) {
    params->BITSEED = malloc(sizeof(*params->BITSEED));
    INITBIT_SEQ(params->BITSEED);
  }
}

#define INIT_RANDSEED_SEQ CONCATENATE(ADDSUFFIXTO(init_randseed), _seq)
static inline void INIT_RANDSEED_SEQ (FPPARAMS *params) {
  if (params->RANDSEED == NULL) {
    params->RANDSEED = malloc(sizeof(*params->RANDSEED));
    INITRAND_SEQ(params->RANDSEED);
  }
}

/* Function to validate floating-point options passed to MAINFUNNAME. */
#define VALIDATE_INPUT CONCATENATE(ADDSUFFIXTO(MAINFUNNAME),_validate_optstruct)
static inline int VALIDATE_INPUT(const optstruct *fpopts) {

  int retval;

  // Set retval to -1 if format is not valid.
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

  //Return 2 if precision is invalid (either nonpositive or too large).
  if (fpopts->precision > DEFPREC || fpopts->precision <= 0)
    return 2;

  // Set retval to -2 if double rounding may occur.
  if (fpopts->precision > (fpopts->round<=1 ? floor((DEFPREC-2)/2) : DEFPREC-1))
    retval = -2;

  // Return 3 if emax is invalid (either nonpositive or too large).
  if (fpopts->emax > DEFEMAX || fpopts->emax <= 0)
    return 3;

  // Set retval to -4 if rounding mode is set to no rounding.
  if (fpopts->round == CPFLOAT_NO_RND)
    retval = -4;

  // Return 5 if p is required but is not a valid probability.
  if (fpopts->flip == CPFLOAT_SOFTERR && (fpopts->p > 1 || fpopts->p < 0))
    return 5;

  // Return 0 or warning value.
  return retval;
}

/* Compute floating-point parameters required by the rounding functions. */
#define COMPUTE_GLOBAL_PARAMS ADDSUFFIXTO(compute_golbal_params)
static inline FPPARAMS COMPUTE_GLOBAL_PARAMS(const optstruct *fpopts,
                                             int *retval) {

  // Actual precision and exponent range.
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

  // Derived floating point parameters.
  const int emin = 1-emax;
  const FPTYPE xmin = ldexp(1., emin);  // Smallest pos. normal.
  FPTYPE ftzthreshold;
  if (fpopts->subnormal == CPFLOAT_SUBN_USE)
    ftzthreshold = ldexp(1., emin-precision+1); // Smallest pos. subnormal.
  else // (fpopts->subnormal == CPFLOAT_SUBN_RND)
    ftzthreshold = xmin;
  const FPTYPE xmax = ldexp(1., emax) * (2-ldexp(1., 1-precision));
  const FPTYPE xbnd = ldexp(1., emax) * (2-ldexp(1., -precision));

  // Bitmasks.
  const INTTYPE leadmask = FULLMASK << (DEFPREC-precision); // Bits to keep
  const INTTYPE trailmask = leadmask ^ FULLMASK; // Bits to discard.

  FPPARAMS params = {precision, emax, emin, fpopts->subnormal,
                     ftzthreshold, xmin, xmax, xbnd,
                     leadmask, trailmask, fpopts->BITSEED, fpopts->RANDSEED};

  return params;
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

/* Round-to-nearest with ties-to-away. */
#define RN_TIES_TO_AWAY_SCALAR_SAME_EXP(x, y, p)                 \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {    \
    if (ABS(y) >= p->xmin/2)                                     \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmin));                \
    else                                                         \
      *(x) = FPOF(SIGN(y));                                      \
  } else                                                         \
    *(x) = FPOF((INTOF(y) +                                      \
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
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                          \
    *(x) = FPOF((INTOF(y) +                                                    \
                 (INTCONST(1) << (DEFPREC-1-lp.precision))) & lp.leadmask);    \
  }

/* Round-to-nearest with ties-to-zero. */
#define RN_TIES_TO_ZERO_SCALAR_SAME_EXP(x, y, p)                 \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {    \
    if (ABS(y) > p->xmin/2)                                      \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmin));                \
    else                                                         \
      *(x) = FPOF(SIGN(y));                                      \
  } else                                                         \
    *(x) = FPOF((INTOF(y) + (p->trailmask>>1)) & p->leadmask);

#define RN_TIES_TO_ZERO_SCALAR_OTHER_EXP(x, y, p, lp)             \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                 \
    if (ABS(y) > p->ftzthreshold/2)                               \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));         \
    else                                                          \
      *(x) = FPOF(SIGN(y));                                       \
  } else if (ABS(y) > p->xbnd) { /* Overflow */                   \
    *(x) = FPOF(SIGN(y) | INTOFCONST(INFINITY));                  \
  } else {                                                        \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                             \
    *(x) = FPOF((INTOF(y) + (lp.trailmask>>1)) & lp.leadmask);    \
  }

/* Round-to-nearest with ties-to-even. */
#define RN_TIES_TO_EVEN_SCALAR_SAME_EXP(x, y, p)                          \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {             \
    if (ABS(y) >= p->xmin/2)                                              \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmin));                         \
    else                                                                  \
      *(x) = FPOF(SIGN(y));                                               \
  } else {                                                                \
    INTTYPE LSB = (INTOF(y) >> (DEFPREC-p->precision)) & INTCONST(1);     \
    *(x) = FPOF((INTOF(y) + (p->trailmask >> 1) + LSB) & p->leadmask);    \
  }

#define RN_TIES_TO_EVEN_SCALAR_OTHER_EXP(x, y, p, lp)                          \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                              \
    if (ABS(y) < p->ftzthreshold/2                                             \
        || (ABS(y) == p->ftzthreshold/2                                        \
            && p->subnormal == CPFLOAT_SUBN_USE))                              \
      *(x) = FPOF(SIGN(y));                                                    \
    else                                                                       \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));                      \
  } else if (ABS(y) >= p->xbnd) { /* Overflow */                               \
    *(x) = FPOF(SIGN(y) | INTOFCONST(INFINITY));                               \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                          \
    INTTYPE LSB = ((INTOF(y) >> (DEFPREC-lp.precision)) & INTCONST(1))         \
      | (lp.precision == 1 && DEFEMAX != p->emax); /* Hidden bit is one. */    \
    *(x) = FPOF((INTOF(y) + ((lp.trailmask >> 1) + LSB)) & lp.leadmask);       \
  }

/* Round-to-plus-infinity (also known as round-up). */
#define RD_TWD_PINF_SCALAR_SAME_EXP(x, y, p)                     \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {    \
    if(*(y) > 0)                                                 \
      *(x) = p->xmin;                                            \
    else                                                         \
      *(x) = FPOF(SIGN(y));                                      \
  } else {                                                       \
    if (*(y) > 0) /* Add ulp if x is positive finite. */         \
      *(x) = FPOF((INTOF(y) + p->trailmask) & p->leadmask);      \
    else                                                         \
      *(x) = FPOF(INTOF(y) & p->leadmask);                       \
  }

#define RD_TWD_PINF_SCALAR_OTHER_EXP(x, y, p, lp)              \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */              \
    *(x) = *(y) > 0 ? p->ftzthreshold : 0;                     \
  } else if (ABS(y) > p->xmax) { /* Overflow */                \
    if (*(y) > p->xmax)                                        \
      *(x) = INFINITY;                                         \
    else if (*(y) < -p->xmax && *(y) != -INFINITY)             \
      *(x) = -p->xmax;                                         \
    else /* *(y) == -INFINITY */                               \
      *(x) = -INFINITY;                                        \
  } else {                                                     \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                          \
    if (SIGN(y) == 0) /* Add ulp if x is positive.  */         \
      *(x) = FPOF((INTOF(y) + lp.trailmask) & lp.leadmask);    \
    else                                                       \
      *(x) = FPOF(INTOF(y) & lp.leadmask);                     \
  }

/* Round-to-minus-infinity (also known as round-down). */
#define RD_TWD_NINF_SCALAR_SAME_EXP(x, y, p)                     \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin) {    \
    if(*(y) < 0)                                                 \
      *(x) = -p->xmin;                                           \
    else                                                         \
      *(x) = FPOF(SIGN(y));                                      \
  } else {                                                       \
    if (*(y) < 0) /* Subtract ulp if x is positive finite. */    \
      *(x) = FPOF((INTOF(y) + p->trailmask) & p->leadmask);      \
    else                                                         \
      *(x) = FPOF(INTOF(y) & p->leadmask);                       \
  }

#define RD_TWD_NINF_SCALAR_OTHER_EXP(x, y, p, lp)              \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */              \
    *(x) = *(y) >= 0 ? 0 : -p->ftzthreshold;                   \
  } else if (ABS(y) > p->xmax) { /* Overflow */                \
    if (*(y) < -p->xmax)                                       \
      *(x) = -INFINITY;                                        \
    else if (*(y) > p->xmax && *(y) != INFINITY)               \
      *(x) = p->xmax;                                          \
    else /* *(y) == INFINITY */                                \
      *(x) = INFINITY;                                         \
  } else {                                                     \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                          \
    if (SIGN(y)) /* Subtract ulp if x is positive. */          \
      *(x) = FPOF((INTOF(y) + lp.trailmask) & lp.leadmask);    \
    else                                                       \
      *(x) = FPOF(INTOF(y) & lp.leadmask);                     \
  }

/* Round-to-zero (also known as truncation). */
#define RD_TWD_ZERO_SCALAR_SAME_EXP(x, y, p)                   \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(y) < p->xmin)    \
    *(x) = FPOF(SIGN(y));                                      \
  else                                                         \
    *(x) = FPOF(INTOF(y) & p->leadmask);

#define RD_TWD_ZERO_SCALAR_OTHER_EXP(x, y, p, lp)                       \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                       \
    *(x) = FPOF(SIGN(y));                                               \
  } else if (ABS(y) > p->xmax && ABS(y) != INFINITY) { /* Overflow */   \
    *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmax));                         \
  } else {                                                              \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                   \
    *(x) = FPOF(INTOF(y) & lp.leadmask);                                \
  }

/* Stochastic rounding with proportional probabilities. */
#define RS_PROP_SCALAR(x, y, p, lp)                                   \
  if (ABS(y) < p->ftzthreshold) { /* Underflow */                     \
    int expx = ((EXPMASK & INTOF(y)) >> (DEFPREC - 1)) - DEFEMAX;     \
    INTTYPE trailfrac = (INTOF(y) & (FULLMASK >> NLEADBITS));         \
    if (expx == -DEFEMAX) { /* No implicit bit (x is subnormal). */   \
      expx += 1;                                                      \
    } else { /* Make implicit bit explicit (x is normal). */          \
      trailfrac = (INTOF(y) & (FULLMASK >> NLEADBITS))                \
        | (INTCONST(1) << (DEFPREC-1));                               \
    }                                                                 \
    int localemin = p->subnormal == CPFLOAT_SUBN_USE ?                \
      p->emin - (int)p->precision + 1: p->emin;                       \
    int expdiff = localemin - expx;                                   \
    INTTYPE rnd = GENRAND(p->RANDSEED) >> 1;                          \
    /* Shift fraction of *(y) left or right as needed. */             \
    if (expdiff <= NLEADBITS - 1)                                     \
      trailfrac <<= NLEADBITS - 1 - expdiff;                          \
    else {                                                            \
      if (expdiff - (NLEADBITS - 1) >= NBITS)                         \
        trailfrac = 0;                                                \
      else                                                            \
        trailfrac >>= (expdiff - (NLEADBITS - 1));                    \
    }                                                                 \
    if (trailfrac > rnd) {                                            \
      *(x) = FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));             \
    } else {                                                          \
      *(x) = FPOF(SIGN(y));                                           \
      continue;                                                       \
    }                                                                 \
  } else if (ABS(y) < p->xbnd) { /* Rounding possibly required. */    \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                 \
    INTTYPE rndbuf = GENRAND(p->RANDSEED) & lp.trailmask;             \
    if (ABS(y) > p->xmax) {                                           \
      rndbuf = rndbuf >> 1;                                           \
      lp.leadmask = (lp.leadmask >> 1) | SIGNMASK;                    \
    }                                                                 \
    *(x) = FPOF((INTOF(y) + rndbuf) & lp.leadmask);                   \
  } else {                                                            \
    *(x) = *(y);                                                      \
  }                                                                   \
  if (ABS(x) >= p->xbnd) /* Overflow */                               \
    *(x) = FPOF(SIGN(y) | INTOFCONST(INFINITY));

/* Stochastic rounding with equal probabilities. */
#define RS_EQUI_SCALAR(x, y, p, lp)                                           \
  UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                           \
  if (ABS(y) < p->ftzthreshold && *(y) != 0) { /* Underflow */                \
    randombit = GENBIT(p->BITSEED);                                           \
    *(x) = FPOF( SIGN(y) | INTOFCONST(randombit ? p->ftzthreshold : 0));      \
  } else if (ABS(y) > p->xmax && ABS(y) != INFINITY) { /* Overflow */         \
    randombit = GENBIT(p->BITSEED);                                           \
    *(x) = FPOF( SIGN(y) | INTOFCONST(randombit ? INFINITY : p->xmax));       \
  } else if ((INTOF(y) & lp.trailmask)) { /* Not exactly representable. */    \
    randombit = GENBIT(p->BITSEED);                                           \
    *(x) = FPOF(INTOF(y) & lp.leadmask);                                      \
    if (randombit)                                                            \
      *(x) = FPOF(INTOF(x) + (INTCONST(1) << (DEFPREC-lp.precision)));        \
  } else /* *(y) exactly representable, no rounding necessary. */             \
    *(x) = *(y);

/* Round-to-odd. */
#define RO_SCALAR(x, y, p, lp)                                                 \
  if (ABS(y) < p->ftzthreshold && *(y) != 0) { /* Underflow */                 \
    *(x) =  FPOF(SIGN(y) | INTOFCONST(p->ftzthreshold));                       \
  } else if (ABS(y) > p->xmax && ABS(y) != INFINITY) { /* Overflow */          \
    *(x) = FPOF(SIGN(y) | INTOFCONST(p->xmax));                                \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                          \
    if ((lp.trailmask & INTOF(y)) /* Not exactly representable. */             \
        && (lp.precision > 1)) /* Set last bit to 1 if stored explicitly. */   \
      *(x) = FPOF((INTOF(y) & lp.leadmask) |                                   \
                  (INTCONST(1) << (DEFPREC-lp.precision)));                    \
    else                                                                       \
      *(x) = FPOF(INTOF(y) & lp.leadmask);                                     \
  }

/******************************
 * Macros for vector rounding *
 ******************************/

#define FOR_STRING(parallel) FOR_STRING_##parallel
#define FOR_STRING_SEQ
#define FOR_STRING_PAR  _Pragma("omp for")

/* Round-to-nearest with ties-to-away. */
#define RN_TIES_TO_AWAY(X, Y, PREPROC, POSTPROC, numelem, p, parallel)    \
  if (p->emax == DEFEMAX) {                                               \
    FOR_STRING(parallel)                                                  \
      for (size_t i=0; i<numelem; i++) {                                  \
        PREPROC                                                           \
        RN_TIES_TO_AWAY_SCALAR_SAME_EXP(X+i, Y+i, p)                      \
        POSTPROC                                                          \
      }                                                                   \
  } else {                                                                \
    FOR_STRING(parallel)                                                  \
      for (size_t i=0; i<numelem; i++) {                                  \
        PREPROC                                                           \
        RN_TIES_TO_AWAY_SCALAR_OTHER_EXP(X+i, Y+i, p, lp)                 \
        POSTPROC                                                          \
      }                                                                   \
  }

/* Round-to-nearest with ties-to-zero. */
#define RN_TIES_TO_ZERO(X, Y, PREPROC, POSTPROC, numelem, p, parallel)    \
  if (p->emax == DEFEMAX) {                                               \
    FOR_STRING(parallel)                                                  \
      for (size_t i=0; i<numelem; i++) {                                  \
        PREPROC                                                           \
        RN_TIES_TO_ZERO_SCALAR_SAME_EXP(X+i, Y+i, p)                      \
        POSTPROC                                                          \
      }                                                                   \
  } else {                                                                \
    FOR_STRING(parallel)                                                  \
      for (size_t i=0; i<numelem; i++) {                                  \
        PREPROC                                                           \
        RN_TIES_TO_ZERO_SCALAR_OTHER_EXP(X+i, Y+i, p, lp)                 \
        POSTPROC                                                          \
      }                                                                   \
  }

/* Round-to-nearest with ties-to-even. */
#define RN_TIES_TO_EVEN(X, Y, PREPROC, POSTPROC, numelem, p, parallel)    \
  if (p->emax == DEFEMAX) {                                               \
    FOR_STRING(parallel)                                                  \
      for (size_t i=0; i<numelem; i++) {                                  \
        PREPROC                                                           \
        RN_TIES_TO_EVEN_SCALAR_SAME_EXP(X+i, Y+i, p)                      \
        POSTPROC                                                          \
      }                                                                   \
  } else {                                                                \
    FOR_STRING(parallel)                                                  \
      for (size_t i=0; i<numelem; i++) {                                  \
        PREPROC                                                           \
        RN_TIES_TO_EVEN_SCALAR_OTHER_EXP(X+i, Y+i, p, lp)                 \
        POSTPROC                                                          \
      }                                                                   \
  }

/* Round-to-plus-infinity (also known as round-up). */
#define RD_TWD_PINF(X, Y, PREPROC, POSTPROC, numelem, p, parallel)    \
  if (p->emax == DEFEMAX) {                                           \
    FOR_STRING(parallel)                                              \
      for (size_t i=0; i<numelem; i++) {                              \
        PREPROC                                                       \
        RD_TWD_PINF_SCALAR_SAME_EXP(X+i, Y+i, p)                      \
        POSTPROC                                                      \
      }                                                               \
  } else {                                                            \
    FOR_STRING(parallel)                                              \
      for (size_t i=0; i<numelem; i++) {                              \
        PREPROC                                                       \
        RD_TWD_PINF_SCALAR_OTHER_EXP(X+i, Y+i, p, lp)                 \
        POSTPROC                                                      \
      }                                                               \
  }

/* Round-to-minus-infinity (also known as round-down). */
#define RD_TWD_NINF(X, Y, PREPROC, POSTPROC, numelem, p, parallel)   \
  if (p->emax == DEFEMAX) {                                          \
    FOR_STRING(parallel)                                             \
      for (size_t i=0; i<numelem; i++) {                             \
        PREPROC                                                      \
        RD_TWD_NINF_SCALAR_SAME_EXP(X+i, Y+i, p)                     \
        POSTPROC                                                     \
      }                                                              \
  } else {                                                           \
    FOR_STRING(parallel)                                             \
      for (size_t i=0; i<numelem; i++) {                             \
        PREPROC                                                      \
        RD_TWD_NINF_SCALAR_OTHER_EXP(X+i, Y+i, p, lp)                \
        POSTPROC                                                     \
      }                                                              \
  }

/* Round-to-zero (also known as truncation). */
#define RD_TWD_ZERO(X, Y, PREPROC, POSTPROC, numelem, p, parallel)    \
  if (p->emax == DEFEMAX) {                                           \
    FOR_STRING(parallel)                                              \
      for (size_t i=0; i<numelem; i++) {                              \
        PREPROC                                                       \
        RD_TWD_ZERO_SCALAR_SAME_EXP (X+i, Y+i, p)                     \
        POSTPROC                                                      \
      }                                                               \
  } else {                                                            \
    FOR_STRING(parallel)                                              \
      for (size_t i=0; i<numelem; i++) {                              \
        PREPROC                                                       \
        RD_TWD_ZERO_SCALAR_OTHER_EXP (X+i, Y+i, p, lp)                \
        POSTPROC                                                      \
      }                                                               \
  }

/* Stochastic rounding with proportional probabilities. */
#define RS_PROP(X, Y, PREPROC, POSTPROC, numelem, p, parallel, INIT_RANDSEED)  \
  INIT_RANDSEED(p);                                                            \
  FOR_STRING(parallel)                                                         \
    for (size_t i=0; i<numelem; i++) {                                         \
      PREPROC                                                                  \
      RS_PROP_SCALAR(X+i, Y+i, p, lp)                                          \
      POSTPROC                                                                 \
    }

/* Stochastic rounding with equal probabilities. */
#define RS_EQUI(X, Y, PREPROC, POSTPROC, numelem, p, parallel, INIT_BITSEED)   \
  INIT_BITSEED(p);                                                             \
  BITTYPE randombit;                                                           \
  FOR_STRING(parallel)                                                         \
    for (size_t i=0; i<numelem; i++) {                                         \
      PREPROC                                                                  \
      RS_EQUI_SCALAR(X+i, Y+i, p, lp)                                          \
      POSTPROC                                                                 \
    }

/* Round-to-odd. */
#define RO(X, Y, PREPROC, POSTPROC, numelem, p, parallel)    \
  FOR_STRING(parallel)                                       \
    for (size_t i=0; i<numelem; i++) {                       \
      PREPROC                                                \
      RO_SCALAR(X+i, Y+i, p, lp)                             \
      POSTPROC                                               \
    }

/********************************************************
 * Macros that define rounding or arithmetic functions. *
 ********************************************************/

#define GENERATE_FUN_NAME(FUNNAME)                              \
  ADDSUFFIXTO(CONCATENATE(CONCATENATE(MAINFUNNAME,_), FUNNAME))
#define GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX)          \
  CONCATENATE(GENERATE_FUN_NAME(FUNNAME), PARALLEL_SUFFIX)

#define PARALLEL_STRING(parallel) PARALLEL_STRING_##parallel
#define PARALLEL_STRING_SEQ
#define PARALLEL_STRING_PAR _Pragma("omp parallel shared(X, A, fpopts)")

/* Macros that generates CPFloat functions. They will have name:
 *     <MAINFUNNAME>_<FUNNAME><TYPESUFFIX><PARALLEL_SUFFIX>
 * The argument list will include FPTYPE *X and the FINPUT string.
 * PREPROC and POSTPROC should be the operations to be performed before
 * and after rounding, respectively.
 */
#define GENERATE_INTERFACE_SUBFUNCTIONS(FUNNAME,                             \
                                        FINPUT,                              \
                                        PREPROC,                             \
                                        POSTPROC,                            \
                                        PARALLEL_TYPE,                       \
                                        PARALLEL_SUFFIX,                     \
                                        INIT_RAND_STRING,                    \
                                        INIT_BIT_STRING)                     \
  static inline int                                                          \
  GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX)                             \
  (FPTYPE *X, FPTYPE *Y, DEPARENTHESIZE_MAYBE(FINPUT),                       \
   const size_t numelem, const optstruct *fpopts) {                          \
    int retval = 0;                                                          \
    FPPARAMS params = COMPUTE_GLOBAL_PARAMS(fpopts, &retval);                \
    FPPARAMS *paramsptr = &params;                                           \
    LOCPARAMS lp;                                                            \
    PARALLEL_STRING(PARALLEL_TYPE)                                           \
      {                                                                      \
        switch (fpopts->round) {                                             \
        case CPFLOAT_RND_NA: /* round-to-nearest with ties-to-away  */       \
          RN_TIES_TO_AWAY(X, Y, PREPROC, POSTPROC, numelem, paramsptr,       \
                          PARALLEL_TYPE);                                    \
          break;                                                             \
        case CPFLOAT_RND_NZ: /* round-to-nearest with ties-to-zero */        \
          RN_TIES_TO_ZERO(X, Y, PREPROC, POSTPROC, numelem, paramsptr,       \
                          PARALLEL_TYPE);                                    \
          break;                                                             \
        case CPFLOAT_RND_NE: /* round-to-nearest with ties-to-even */        \
          RN_TIES_TO_EVEN(X, Y, PREPROC, POSTPROC, numelem, paramsptr,       \
                          PARALLEL_TYPE);                                    \
          break;                                                             \
        case CPFLOAT_RND_TP: /* round-toward-positive */                     \
          RD_TWD_PINF(X, Y, PREPROC, POSTPROC, numelem, paramsptr,           \
                      PARALLEL_TYPE);                                        \
          break;                                                             \
        case CPFLOAT_RND_TN: /* round-toward-negative */                     \
          RD_TWD_NINF(X, Y, PREPROC, POSTPROC, numelem, paramsptr,           \
                      PARALLEL_TYPE);                                        \
          break;                                                             \
        case CPFLOAT_RND_TZ: /* round-toward-zero */                         \
          RD_TWD_ZERO(X, Y, PREPROC, POSTPROC, numelem, paramsptr,           \
                      PARALLEL_TYPE);                                        \
          break;                                                             \
        case CPFLOAT_RND_SP: /* stochastic rounding (proportional) */        \
          RS_PROP(X, Y, PREPROC, POSTPROC, numelem, paramsptr,               \
                  PARALLEL_TYPE, INIT_RAND_STRING);                          \
          break;                                                             \
        case CPFLOAT_RND_SE: /* stochastic rounding (equal) */               \
          RS_EQUI(X, Y, PREPROC, POSTPROC, numelem, paramsptr,               \
                  PARALLEL_TYPE, INIT_BIT_STRING);                           \
          break;                                                             \
        case CPFLOAT_RND_OD: /* round-to-odd */                              \
          RO(X, Y, PREPROC, POSTPROC, numelem, paramsptr, PARALLEL_TYPE);    \
          break;                                                             \
        default: /* No rounding if unknown mode specified. */                \
          FOR_STRING(PARALLEL_TYPE)                                          \
            for (size_t i=0; i<numelem; i++) {                               \
              X[i] = A[i];                                                   \
            }                                                                \
          break;                                                             \
        }                                                                    \
      }                                                                      \
    /* Introduce bit flips. */                                               \
    if (fpopts->flip == CPFLOAT_SOFTERR) {                                   \
      FOR_STRING(PARALLEL_TYPE)                                              \
        for (size_t i=0; i<numelem; i++) {                                   \
          if (rand() / (FPTYPE)RAND_MAX < fpopts->p) {                       \
            X[i] = FPOF(INTOF(X+i) ^                                         \
                        (INTCONST(1) << rand() % (DEFPREC - 1)));            \
          }                                                                  \
        }                                                                    \
    }                                                                        \
    return retval;                                                           \
  }

#define PARALLEL_SUFFIX_SEQ _sequential
#define INIT_RAND_STRING_SEQ(p) INIT_RANDSEED_SEQ(p);
#define INIT_BIT_STRING_SEQ(p) INIT_BITSEED_SEQ(p);

#define GENERATE_INTERFACE_SUBFUNCTIONS_SEQ(FUNNAME, FINPUT,      \
                                            PREPROC, POSTPROC)    \
  GENERATE_INTERFACE_SUBFUNCTIONS(FUNNAME,                        \
                                  FINPUT,                         \
                                  PREPROC,                        \
                                  POSTPROC,                       \
                                  SEQ,                            \
                                  PARALLEL_SUFFIX_SEQ,            \
                                  INIT_RAND_STRING_SEQ,           \
                                  INIT_BIT_STRING_SEQ)

#ifdef _OPENMP
#define PARALLEL_SUFFIX_PAR _parallel
#define INIT_RAND_STRING_PAR INIT_RANDSEED_PAR
#define INIT_BIT_STRING_PAR INIT_BITSEED_PAR

#define GENERATE_INTERFACE_SUBFUNCTIONS_PAR(FUNNAME, FINPUT,      \
                                            PREPROC, POSTPROC)    \
  GENERATE_INTERFACE_SUBFUNCTIONS(FUNNAME,                        \
                                  FINPUT,                         \
                                  PREPROC,                        \
                                  POSTPROC,                       \
                                  PAR,                            \
                                  PARALLEL_SUFFIX_PAR,            \
                                  INIT_RAND_STRING_PAR,           \
                                  INIT_BIT_STRING_PAR)

// Partial results before rounding are store in the vector X.
#define GENERATE_INTERFACE(FUNNAME, Y, FINPUT, FCALL, PREPROC, POSTPROC)     \
  GENERATE_INTERFACE_SUBFUNCTIONS_SEQ(FUNNAME, FINPUT, PREPROC, POSTPROC)    \
  GENERATE_INTERFACE_SUBFUNCTIONS_PAR(FUNNAME, FINPUT, PREPROC, POSTPROC)    \
  static inline int                                                          \
  GENERATE_FUN_NAME(FUNNAME)(FPTYPE *X,                                      \
                             DEPARENTHESIZE_MAYBE(FINPUT),                   \
                             const size_t numelem,                           \
                             const optstruct *fpopts) {                      \
    if (numelem < CONCATENATE(OPENMP_THRESHOLD_, FPTYPE))                    \
      return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_SEQ)              \
        (X, Y, DEPARENTHESIZE_MAYBE(FCALL), numelem, fpopts);                \
    else                                                                     \
      return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_PAR)              \
        (X, Y, DEPARENTHESIZE_MAYBE(FCALL), numelem, fpopts);                \
  }
#else /* #ifdef _OPENMP */
#define GENERATE_INTERFACE(FUNNAME, Y, FINPUT, FCALL, PREPROC, POSTPROC)     \
  GENERATE_INTERFACE_SUBFUNCTIONS_SEQ(FUNNAME, FINPUT, PREPROC, POSTPROC)    \
  static inline int                                                          \
       GENERATE_FUN_NAME(FUNNAME)(FPTYPE *X,                                 \
                             DEPARENTHESIZE_MAYBE(FINPUT),                   \
                             const size_t numelem,                           \
                             const optstruct *fpopts) {                      \
    return GENERATE_SUBFUN_NAME(FUNNAME, PARALLEL_SUFFIX_SEQ)                \
      (X, Y, DEPARENTHESIZE_MAYBE(FCALL), numelem, fpopts);                  \
  }
#endif /* #ifdef _OPENMP */

/**************************************
 * Actual instantiation of functions. *
 **************************************/

/* Only rounding. */
GENERATE_INTERFACE(round, (FPTYPE *)A, const FPTYPE* A, A, NOOP, NOOP)
static inline
int ADDSUFFIXTO(cpfloat)(FPTYPE *X,
                         const FPTYPE *A,
                         const size_t numelem,
                         const optstruct *fpopts) {
  return GENERATE_FUN_NAME(round)(X, A, numelem, fpopts);
}


/* Undefine macros from this file. */

#undef NOOP

#undef DO_REMOVE_PARENTHESES
#undef EXPAND_INTERIOR
#undef REMOVE_INTERIOR_
#undef REMOVE_PARENTHESES
#undef DEPARENTHESIZE_MAYBE

#undef CONCATENATE_INNER
#undef CONCATENATE

#undef BITSEED
#undef BITSEEDTYPE
#undef BITTYPE
#undef INITBIT_SEQ
#undef INITBIT_PAR
#undef GENBIT

#undef MAINFUNNAME

#undef ADDSUFFIXTO
#undef FPUNION
#undef RANDSEED
#undef RANDSEEDTYPE
#undef FPPARAMS
#undef LOCPARAMS
#undef INTCONST
#undef INTOF
#undef INTOFCONST
#undef FPOF
#undef SIGN
#undef ABS

#undef INIT_RANDSEED_PAR
#undef INIT_BITSEED_PAR
#undef INIT_BITSEED_SEQ
#undef INIT_RANDSEED_SEQ

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

#undef INITRAND_SEQ
#undef INITRAND_PAR
#undef GENRAND

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
