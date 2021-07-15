/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * This file is part of CPFloat.
 *
 * The macros and functions defined here require a number of macros be defined
 * before this header is included. This file should never be included directly
 * in production code. In order to use CPFloat, it suffices to include
 * cpfloat_binary32.h and cpfloat_binary64.h. These two headers include this
 * file multiple times handling the definition of macros correctly.
 */

/**************************************
 **************************************
 * MACROS DEFINED ONLY ONCE PER TYPE. *
 **************************************
 **************************************/
#ifndef CONCATENATE_INNER
#define CONCATENATE_INNER(arg1, arg2) arg1 ## arg2
#define CONCATENATE(arg1, arg2) CONCATENATE_INNER(arg1, arg2)
#define ADDSUFFIXTO(x) CONCATENATE(x, FUNSUFFIX)
#define MAINFUNNAME cpfloat

/* Struct to reinterpret floating-point values as unsigned integers. */
#define FPUNION ADDSUFFIXTO(fpint)
typedef union {
  INTTYPE intval;
  FPTYPE fpval;
} FPUNION;

/* Types for internal state of pseudo-random number generator. */
#define BITSEED bitseed
#define BITSEEDTYPE cpfloat_bitseed_t
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

/* Functions to initialize bit pseudo-random generators and generate bits. */
/*
 * NOTE: The following block does not depend on the type of the array (float or
 * double), thus it could in principle be moved out of cpfloat_template.h.
 */
#ifdef PCG_VARIANTS_H_INCLUDED
#define BITTYPE unsigned int
#define INITBIT_SINGLE(seed)                    \
  pcg32_srandom_r(seed,                         \
                  time(NULL),                   \
                  (intptr_t)seed)
#define INITBIT_MULTI(seed)                                     \
  pcg32_srandom_r(seed,                                         \
                  omp_get_thread_num() * 13254 + time(NULL),    \
                  (intptr_t)seed)
#define GENBIT(seed) (pcg32_random_r(seed) & (1U << 31))
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#define BITTYPE unsigned int
#define INITBIT_SINGLE(seed) *seed = time(NULL)
#define INITBIT_MULTI(seed) *seed = omp_get_thread_num() * 13254 + time(NULL)
#define GENBIT(seed) (rand_r(seed) & (1U << 30))
#endif /* #ifdef PCG_VARIANTS_H_INCLUDED */

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
#define RN_TIES_TO_AWAY_SCALAR_SAME_EXP(x, a, p)                \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(a) < p->xmin) {    \
    if (ABS(a) >= p->xmin/2)                                     \
      *(x) = FPOF(SIGN(a) | INTOFCONST(p->xmin));              \
    else                                                           \
      *(x) = FPOF(SIGN(a));                                    \
  } else                                                           \
    *(x) = FPOF((INTOF(a) +                                    \
                 (INTCONST(1) << (DEFPREC-1-p->precision))) & p->leadmask);

#define RN_TIES_TO_AWAY_SCALAR_OTHER_EXP(x, a, p, lp)                       \
  if (ABS(a) < p->ftzthreshold) { /* Underflow */                            \
   if (ABS(a) < p->ftzthreshold/2)                                           \
     *(x) = FPOF(SIGN(a));                                                 \
   else                                                                        \
     *(x) = FPOF(SIGN(a) | INTOFCONST(p->ftzthreshold));                   \
  } else if (ABS(a) >= p->xbnd) { /* Overflow */                             \
    *(x) = FPOF(SIGN(a) | INTOFCONST(INFINITY));                           \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                          \
    *(x) = FPOF((INTOF(a) +                                                \
                 (INTCONST(1) << (DEFPREC-1-lp.precision))) & lp.leadmask);    \
  }

/* Round-to-nearest with ties-to-zero. */
#define RN_TIES_TO_ZERO_SCALAR_SAME_EXP(x, a, p)                \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(a) < p->xmin) {    \
    if (ABS(a) > p->xmin/2)                                      \
      *(x) = FPOF(SIGN(a) | INTOFCONST(p->xmin));              \
    else                                                           \
      *(x) = FPOF(SIGN(a));                                    \
  } else                                                           \
    *(x) = FPOF((INTOF(a) + (p->trailmask>>1)) & p->leadmask);

#define RN_TIES_TO_ZERO_SCALAR_OTHER_EXP(x, a, p, lp)              \
  if (ABS(a) < p->ftzthreshold) { /* Underflow */                   \
    if (ABS(a) > p->ftzthreshold/2)                                 \
      *(x) = FPOF(SIGN(a) | INTOFCONST(p->ftzthreshold));         \
    else                                                              \
      *(x) = FPOF(SIGN(a));                                       \
  } else if (ABS(a) > p->xbnd) { /* Overflow */                     \
    *(x) = FPOF(SIGN(a) | INTOFCONST(INFINITY));                  \
  } else {                                                            \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                 \
    *(x) = FPOF((INTOF(a) + (lp.trailmask>>1)) & lp.leadmask);    \
  }

/* Round-to-nearest with ties-to-even. */
#define RN_TIES_TO_EVEN_SCALAR_SAME_EXP(x, a, p)                           \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(a) < p->xmin) {               \
    if (ABS(a) >= p->xmin/2)                                                \
          *(x) = FPOF(SIGN(a) | INTOFCONST(p->xmin));                     \
    else                                                                      \
          *(x) = FPOF(SIGN(a));                                           \
  } else {                                                                    \
    INTTYPE LSB = (INTOF(a) >> (DEFPREC-p->precision)) & INTCONST(1);       \
    *(x) = FPOF((INTOF(a) + (p->trailmask >> 1) + LSB) & p->leadmask);    \
  }

#define RN_TIES_TO_EVEN_SCALAR_OTHER_EXP(x, a, p, lp)                       \
  if (ABS(a) < p->ftzthreshold) { /* Underflow */                            \
    if (ABS(a) < p->ftzthreshold/2                                           \
        || (ABS(a) == p->ftzthreshold/2                                      \
            && p->subnormal == CPFLOAT_SUBN_USE))                              \
      *(x) = FPOF(SIGN(a));                                                \
    else                                                                       \
      *(x) = FPOF(SIGN(a) | INTOFCONST(p->ftzthreshold));                  \
  } else if (ABS(a) >= p->xbnd) { /* Overflow */                             \
    *(x) = FPOF(SIGN(a) | INTOFCONST(INFINITY));                           \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                          \
    INTTYPE LSB = ((INTOF(a) >> (DEFPREC-lp.precision)) & INTCONST(1))       \
      | (lp.precision == 1 && DEFEMAX != p->emax); /* Hidden bit is one. */    \
    *(x) = FPOF((INTOF(a) + ((lp.trailmask >> 1) + LSB)) & lp.leadmask);   \
  }

/* Round-to-plus-infinity (also known as round-up). */
#define RD_TWD_PINF_SCALAR_SAME_EXP(x, a, p)                    \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(a) < p->xmin) {    \
    if(*(a) > 0)                                                 \
      *(x) = p->xmin;                                            \
    else                                                           \
      *(x) = FPOF(SIGN(a));                                    \
  } else {                                                         \
    if (*(a) > 0) /* Add ulp if x is positive finite. */         \
      *(x) = FPOF((INTOF(a) + p->trailmask) & p->leadmask);    \
    else                                                           \
      *(x) = FPOF(INTOF(a) & p->leadmask);                     \
  }

#define RD_TWD_PINF_SCALAR_OTHER_EXP(x, a, p, lp)               \
  if (ABS(a) < p->ftzthreshold) { /* Underflow */                \
    *(x) = *(a) > 0 ? p->ftzthreshold : 0;                     \
  } else if (ABS(a) > p->xmax) { /* Overflow */                  \
    if (*(a) > p->xmax)                                          \
      *(x) = INFINITY;                                           \
    else if (*(a) < -p->xmax && *(a) != -INFINITY)             \
      *(x) = -p->xmax;                                           \
    else /* *(a) == -INFINITY */                                 \
      *(x) = -INFINITY;                                          \
  } else {                                                         \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                              \
    if (SIGN(a) == 0) /* Add ulp if x is positive.  */           \
      *(x) = FPOF((INTOF(a) + lp.trailmask) & lp.leadmask);    \
    else                                                           \
      *(x) = FPOF(INTOF(a) & lp.leadmask);                     \
  }

/* Round-to-minus-infinity (also known as round-down). */
#define RD_TWD_NINF_SCALAR_SAME_EXP(x, a, p)                    \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(a) < p->xmin) {    \
    if(*(a) < 0)                                                 \
      *(x) = -p->xmin;                                           \
    else                                                           \
      *(x) = FPOF(SIGN(a));                                    \
  } else {                                                         \
    if (*(a) < 0) /* Subtract ulp if x is positive finite. */    \
      *(x) = FPOF((INTOF(a) + p->trailmask) & p->leadmask);    \
    else                                                           \
      *(x) = FPOF(INTOF(a) & p->leadmask);                     \
  }

#define RD_TWD_NINF_SCALAR_OTHER_EXP(x, a, p, lp)               \
  if (ABS(a) < p->ftzthreshold) { /* Underflow */                \
    *(x) = *(a) >= 0 ? 0 : -p->ftzthreshold;                   \
  } else if (ABS(a) > p->xmax) { /* Overflow */                  \
    if (*(a) < -p->xmax)                                         \
      *(x) = -INFINITY;                                          \
    else if (*(a) > p->xmax && *(a) != INFINITY)               \
      *(x) = p->xmax;                                            \
    else /* *(a) == INFINITY */                                  \
      *(x) = INFINITY;                                           \
  } else {                                                         \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                              \
    if (SIGN(a)) /* Subtract ulp if x is positive. */            \
      *(x) = FPOF((INTOF(a) + lp.trailmask) & lp.leadmask);    \
    else                                                           \
      *(x) = FPOF(INTOF(a) & lp.leadmask);                     \
  }

/* Round-to-zero (also known as truncation). */
#define RD_TWD_ZERO_SCALAR_SAME_EXP(x, a, p)                  \
  if (p->subnormal != CPFLOAT_SUBN_USE && ABS(a) < p->xmin)    \
    *(x) = FPOF(SIGN(a));                                    \
  else                                                           \
    *(x) = FPOF(INTOF(a) & p->leadmask);

#define RD_TWD_ZERO_SCALAR_OTHER_EXP(x, a, p, lp)                         \
  if (ABS(a) < p->ftzthreshold) { /* Underflow */                          \
    *(x) = FPOF(SIGN(a));                                                \
  } else if (ABS(a) > p->xmax && ABS(a) != INFINITY) { /* Overflow */    \
    *(x) = FPOF(SIGN(a) | INTOFCONST(p->xmax));                          \
  } else {                                                                   \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                        \
    *(x) = FPOF(INTOF(a) & lp.leadmask);                                 \
  }


/* Stochastic rounding with proportional probabilities. */
#define RS_PROP_SCALAR(x, a, p, lp)                                  \
  if (ABS(a) < p->ftzthreshold){ /* Underflow */                      \
    int expx = ((EXPMASK & INTOF(a)) >> (DEFPREC - 1)) - DEFEMAX;     \
    INTTYPE trailfrac = (INTOF(a) & (FULLMASK >> NLEADBITS));         \
    if (expx == -DEFEMAX) { /* No implicit bit (x is subnormal). */     \
      expx += 1;                                                        \
    } else { /* Make implicit bit explicit (x is normal). */            \
      trailfrac = (INTOF(a) & (FULLMASK >> NLEADBITS))                \
        | (INTCONST(1) << (DEFPREC-1));                                 \
    }                                                                   \
    int localemin = p->subnormal == CPFLOAT_SUBN_USE ?                  \
      p->emin - (int)p->precision + 1: p->emin;                         \
    int expdiff = localemin - expx;                                     \
    INTTYPE rnd = GENRAND(p->RANDSEED) >> 1;                            \
    /* Shift fraction of *(a) left or right as needed. */             \
    if (expdiff <= NLEADBITS - 1)                                       \
      trailfrac <<= NLEADBITS - 1 - expdiff;                            \
    else {                                                              \
      if (expdiff - (NLEADBITS - 1) >= NBITS)                           \
        trailfrac = 0;                                                  \
      else                                                              \
        trailfrac >>= (expdiff - (NLEADBITS - 1));                      \
    }                                                                   \
    if (trailfrac > rnd) {                                              \
      *(x) = FPOF(SIGN(a) | INTOFCONST(p->ftzthreshold));           \
    } else {                                                            \
      *(x) = FPOF(SIGN(a));                                         \
      continue;                                                         \
    }                                                                   \
  } else if (ABS(a) < p->xbnd) { /* Rounding possibly required. */    \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                   \
    INTTYPE rndbuf = GENRAND(p->RANDSEED) & lp.trailmask;               \
    if (ABS(a) > p->xmax) {                                           \
      rndbuf = rndbuf >> 1;                                             \
      lp.leadmask = (lp.leadmask >> 1) | SIGNMASK;                      \
    }                                                                   \
    *(x) = FPOF((INTOF(a) + rndbuf) & lp.leadmask);                 \
  } else {                                                              \
    *(x) = *(a);                                                    \
  }                                                                     \
  if (ABS(x) >= p->xbnd) /* Overflow */                               \
    *(x) = FPOF(SIGN(a) | INTOFCONST(INFINITY));

/* Stochastic rounding with equal probabilities. */
#define RS_EQUI_SCALAR(x, a, p, lp)                                         \
  UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                            \
  if (ABS(a) < p->ftzthreshold && *(a) != 0) { /* Underflow */             \
    randombit = GENBIT(p->BITSEED);                                            \
    *(x) = FPOF( SIGN(a) | INTOFCONST(randombit ? p->ftzthreshold : 0));   \
  } else if (ABS(a) > p->xmax && ABS(a) != INFINITY) { /* Overflow */      \
    randombit = GENBIT(p->BITSEED);                                            \
    *(x) = FPOF( SIGN(a) | INTOFCONST(randombit ? INFINITY : p->xmax));    \
  } else if ((INTOF(a) & lp.trailmask)) { /* Not exactly representable. */   \
    randombit = GENBIT(p->BITSEED);                                            \
    *(x) = FPOF(INTOF(a) & lp.leadmask);                                   \
    if (randombit)                                                             \
      *(x) = FPOF(INTOF(x) + (INTCONST(1) << (DEFPREC-lp.precision)));     \
  } else /* *(a) exactly representable, no rounding necessary. */            \
    *(x) = *(a);

/* Round-to-odd. */
#define RO_SCALAR(x, a, p, lp)                                              \
  if (ABS(a) < p->ftzthreshold && *(a) != 0) { /* Underflow */             \
    *(x) =  FPOF(SIGN(a) | INTOFCONST(p->ftzthreshold));                   \
  } else if (ABS(a) > p->xmax && ABS(a) != INFINITY) { /* Overflow */      \
    *(x) = FPOF(SIGN(a) | INTOFCONST(p->xmax));                            \
  } else {                                                                     \
    UPDATE_LOCAL_PARAMS(A+i, p, &lp);                                          \
    if ((lp.trailmask & INTOF(a)) /* Not exactly representable. */           \
        && (lp.precision > 1)) /* Set last bit to 1 if stored explicitly. */   \
      *(x) = FPOF((INTOF(a) & lp.leadmask) |                               \
                  (INTCONST(1) << (DEFPREC-lp.precision)));                    \
    else                                                                       \
      *(x) = FPOF(INTOF(a) & lp.leadmask);                                 \
  }

#endif /* #ifndef CONCATENATE_INNER */




/************************************************
 ************************************************
 * MACROS WITH SEQUENTIAL AND PARALLEL FLAVOUR. *
 ************************************************
 ************************************************/
/*  Rounding functions. */
#ifdef SINGLE_THREADED
#define RN_TIES_TO_AWAY CONCATENATE(ADDSUFFIXTO(rn_tta), _sequential)
#define RN_TIES_TO_ZERO CONCATENATE(ADDSUFFIXTO(rn_ttz), _sequential)
#define RN_TIES_TO_EVEN CONCATENATE(ADDSUFFIXTO(rn_tte), _sequential)

#define RD_TWD_PINF CONCATENATE(ADDSUFFIXTO(rd_pinf), _sequential)
#define RD_TWD_NINF CONCATENATE(ADDSUFFIXTO(rd_ninf), _sequential)
#define RD_TWD_ZERO CONCATENATE(ADDSUFFIXTO(rd_zero), _sequential)

#define RS_PROP CONCATENATE(ADDSUFFIXTO(rs_prop), _sequential)
#define RS_EQUI CONCATENATE(ADDSUFFIXTO(rs_equi), _sequential)

#define RO CONCATENATE(ADDSUFFIXTO(ro), _sequential)

/* Functions to initialize the pseudo-random number generator. */
#define INIT_BITSEED_SINGLE CONCATENATE(ADDSUFFIXTO(init_bitseed), _single)
static inline void INIT_BITSEED_SINGLE (FPPARAMS *params) {
  if (params->BITSEED == NULL) {
    params->BITSEED = malloc(sizeof(*params->BITSEED));
    INITBIT_SINGLE(params->BITSEED);
  }
}
#define INIT_RANDSEED_SINGLE CONCATENATE(ADDSUFFIXTO(init_seed), _single)
static inline void INIT_RANDSEED_SINGLE (FPPARAMS *params) {
  if (params->RANDSEED == NULL) {
    params->RANDSEED = malloc(sizeof(*params->RANDSEED));
    INITRAND_SINGLE(params->RANDSEED);
  }
}
#else /* #ifdef SINGLE_THREADED */
#define RN_TIES_TO_AWAY CONCATENATE(ADDSUFFIXTO(rn_tta), _parallel)
#define RN_TIES_TO_ZERO CONCATENATE(ADDSUFFIXTO(rn_ttz), _parallel)
#define RN_TIES_TO_EVEN CONCATENATE(ADDSUFFIXTO(rn_tte), _parallel)

#define RD_TWD_PINF CONCATENATE(ADDSUFFIXTO(rd_pinf), _parallel)
#define RD_TWD_NINF CONCATENATE(ADDSUFFIXTO(rd_ninf), _parallel)
#define RD_TWD_ZERO CONCATENATE(ADDSUFFIXTO(rd_zero), _parallel)

#define RS_PROP CONCATENATE(ADDSUFFIXTO(rs_prop), _parallel)
#define RS_EQUI CONCATENATE(ADDSUFFIXTO(rs_equi), _parallel)

#define RO CONCATENATE(ADDSUFFIXTO(ro), _parallel)

/* Functions to initialize the pseudo-random number generator. */
#define INIT_BITSEED_MULTI CONCATENATE(ADDSUFFIXTO(init_bitseed), _multi)
static inline void INIT_BITSEED_MULTI (FPPARAMS *params) {
  if (params->BITSEED == NULL) {
    params->BITSEED = malloc(sizeof(*params->BITSEED));
    INITBIT_MULTI(params->BITSEED);
  }
}
#define INIT_RANDSEED_MULTI CONCATENATE(ADDSUFFIXTO(init_seed), _multi)
static inline void INIT_RANDSEED_MULTI (FPPARAMS *params) {
  if (params->RANDSEED == NULL) {
    params->RANDSEED = malloc(sizeof(*params->RANDSEED));
    INITRAND_MULTI(params->RANDSEED);
  }
}
#endif /* #ifdef SINGLE_THREADED */

/* Routine for round-to-nearest with ties-to-away. */
static inline void RN_TIES_TO_AWAY(FPTYPE *X,
                                   const FPTYPE *A,
                                   const size_t numelem,
                                   const FPPARAMS *p) {
  if (p->emax == DEFEMAX) {
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RN_TIES_TO_AWAY_SCALAR_SAME_EXP(X+i, A+i, p)
    }
  } else {
    LOCPARAMS lp;
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RN_TIES_TO_AWAY_SCALAR_OTHER_EXP(X+i, A+i, p, lp)
    }
  }
}

/* Routine for round-to-nearest with ties-to-zero. */
static inline void RN_TIES_TO_ZERO(FPTYPE *X,
                                   const FPTYPE *A,
                                   const size_t numelem,
                                   const FPPARAMS *p) {
  if (p->emax == DEFEMAX) {
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RN_TIES_TO_ZERO_SCALAR_SAME_EXP(X+i, A+i, p)
    }
  } else {
    LOCPARAMS lp;
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RN_TIES_TO_ZERO_SCALAR_OTHER_EXP(X+i, A+i, p, lp)
    }
  }
}

/* Routine for round-to-nearest with ties-to-even. */
static inline void RN_TIES_TO_EVEN(FPTYPE *X,
                                   const FPTYPE *A,
                                   const size_t numelem,
                                   const FPPARAMS *p) {
  if (p->emax == DEFEMAX) {
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RN_TIES_TO_EVEN_SCALAR_SAME_EXP(X+i, A+i, p)
    }
  } else {
    LOCPARAMS lp;
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RN_TIES_TO_EVEN_SCALAR_OTHER_EXP(X+i, A+i, p, lp)
    }
  }
}

/* Routine for round-to-plus-infinity (also known as round-up). */
static inline void RD_TWD_PINF(FPTYPE *X,
                               const FPTYPE *A,
                               const size_t numelem,
                               const FPPARAMS *p) {
  if (p->emax == DEFEMAX) {
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RD_TWD_PINF_SCALAR_SAME_EXP(X+i, A+i, p)
    }
  } else {
    LOCPARAMS lp;
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RD_TWD_PINF_SCALAR_OTHER_EXP(X+i, A+i, p, lp)
    }
  }
}

/* Routine for round-to-minus-infinity (also known as round-down). */
static inline void RD_TWD_NINF(FPTYPE *X,
                               const FPTYPE *A,
                               const size_t numelem,
                               const FPPARAMS *p) {
  if (p->emax == DEFEMAX) {
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RD_TWD_NINF_SCALAR_SAME_EXP(X+i, A+i, p)
    }
  } else {
    LOCPARAMS lp;
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RD_TWD_NINF_SCALAR_OTHER_EXP(X+i, A+i, p, lp)
    }
  }
}

/* Routine for round-to-zero (also known as truncation). */
static inline void RD_TWD_ZERO(FPTYPE *X,
                               const FPTYPE *A,
                               const size_t numelem,
                               const FPPARAMS *p) {
  if (p->emax == DEFEMAX) {
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RD_TWD_ZERO_SCALAR_SAME_EXP (X+i, A+i, p)
    }
  } else {
    LOCPARAMS lp;
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      RD_TWD_ZERO_SCALAR_OTHER_EXP (X+i, A+i, p, lp)
    }
  }
}

/* Routine for stochastic rounding with proportional probabilities. */
static inline void RS_PROP(FPTYPE *X,
                           const FPTYPE *A,
                           const size_t numelem,
                           FPPARAMS *p) {
  LOCPARAMS lp;
  #ifdef USE_OPENMP
  INIT_RANDSEED_MULTI(p);
  #pragma omp for
  #else /* #ifdef USE_OPENMP */
  INIT_RANDSEED_SINGLE(p);
  #endif /* #ifdef USE_OPENMP */
  for (size_t i=0; i<numelem; i++){
    RS_PROP_SCALAR(X+i, A+i, p, lp)
  }
}

/* Routine for stochastic rounding with equal probabilities. */
static inline void RS_EQUI(FPTYPE *X,
                           const FPTYPE *A,
                           const size_t numelem,
                           FPPARAMS *p) {
  LOCPARAMS lp;
  BITTYPE randombit;
  #ifdef USE_OPENMP
  INIT_BITSEED_MULTI(p);
  #pragma omp for
  #else /* #ifdef USE_OPENMP */
  INIT_BITSEED_SINGLE(p);
  #endif /* #ifdef USE_OPENMP */
  for (size_t i=0; i<numelem; i++){
    RS_EQUI_SCALAR(X+i, A+i, p, lp)
  }
}

/* Routine for round-to-odd. */
static inline void RO(FPTYPE *X,
                      const FPTYPE *A,
                      const size_t numelem,
                      const FPPARAMS *p) {
  LOCPARAMS lp;
  #ifdef USE_OPENMP
  #pragma omp for
  #endif /* #ifdef USE_OPENMP */
  for (size_t i=0; i<numelem; i++){
    RO_SCALAR(X+i, A+i, p, lp)
  }
}

/*
 * Macros that define the main rounding functions. Either two or three functions
 * will be generated for each storage format, depending on whether OpenMP is
 * supported or not.
 *
 * 1. When SINGLE_THREADED is defined, the following code generates the
 *    single-threaded function cpfloatf?_sequential.
 *
 * 2. If _OPENMP is defined (i.e., OpenMP support is available), then the
 *    following code generates the multi-threaded function cpfloatf?_parallel.
 *
 * 3. The function cpfloatf? is generated as MAINFUN_COMBO below.
 */
#ifdef SINGLE_THREADED
#define MAINFUN CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _sequential)
#else /* #ifdef SINGLE_THREADED */
#define MAINFUN CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _parallel)
#define USE_OPENMP
#endif /* #ifdef SINGLE_THREADED */
static inline int MAINFUN(FPTYPE *X,
                          const FPTYPE *A,
                          const size_t numelem,
                          const optstruct *fpopts) {

  int retval = 0;
  FPPARAMS p = COMPUTE_GLOBAL_PARAMS(fpopts, &retval);

  #ifdef USE_OPENMP
  #pragma omp parallel shared(X, A, fpopts)
  #endif /* #ifdef USE_OPENMP */
  {
    switch (fpopts->round) {
    case CPFLOAT_RND_NA: // round-to-nearest with ties-to-away
      RN_TIES_TO_AWAY(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_NZ: // round-to-nearest with ties-to-zero
      RN_TIES_TO_ZERO(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_NE: // round-to-nearest with ties-to-even
      RN_TIES_TO_EVEN(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_TP: // round-toward-positive
      RD_TWD_PINF(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_TN: // round-toward-negative
      RD_TWD_NINF(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_TZ: // round-toward-zero
      RD_TWD_ZERO(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_SP: // stochastic rounding with proportional probabilities
      RS_PROP(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_SE: // stochastic rounding with equal probabilities
      RS_EQUI(X, A, numelem, &p);
      break;
    case CPFLOAT_RND_OD: // round-to-odd
      RO(X, A, numelem, &p);
      break;
    default: // No rounding if unknown mode specified.
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (size_t i=0; i<numelem; i++){
        X[i] = A[i];
      }
      break;
    }
  }

  // Bit flips.
  if (fpopts->flip == CPFLOAT_SOFTERR) {
    #ifdef USE_OPENMP
    #pragma omp for
    #endif /* #ifdef USE_OPENMP */
    for (size_t i=0; i<numelem; i++){
      if (rand() / (FPTYPE)RAND_MAX < fpopts->p) {
        X[i] = FPOF(INTOF(X+i) ^ (INTCONST(1) << rand() % (DEFPREC - 1)));
      }
    }
  }

  return retval;
}

#define MAINFUN_COMBO ADDSUFFIXTO(MAINFUNNAME)
#ifdef _OPENMP
#ifndef USE_OPENMP
#define MAINFUN_SINGLE CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _sequential)
#define MAINFUN_MULTI  CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _parallel)
static inline int MAINFUN_COMBO(FPTYPE *X,
                                const FPTYPE *A,
                                const size_t numelem,
                                const optstruct *fpopts) {
  if (numelem < CONCATENATE(OPENMP_THRESHOLD_, FPTYPE))
    return MAINFUN_SINGLE(X, A, numelem, fpopts);
  else
    return MAINFUN_MULTI(X, A, numelem, fpopts);
}
#undef MAINFUN_SINGLE
#undef MAINFUN_MULTI
#undef MAINFUN_COMBO
#endif /*#ifdef USE_OPENMP */
#else /* #ifdef _OPENMP */
static inline int MAINFUN_COMBO(FPTYPE *X,
                                const FPTYPE *A,
                                const size_t numelem,
                                const optstruct *fpopts) {
  return MAINFUN(X, A, numelem, fpopts);
}
#endif /* #ifdef _OPENMP */

#undef RN_TIES_TO_AWAY
#undef RN_TIES_TO_ZERO
#undef RN_TIES_TO_EVEN
#undef RD_TWD_PINF
#undef RD_TWD_NINF
#undef RD_TWD_ZERO
#undef RS_PROP
#undef RS_EQUI
#undef RO
#undef MAINFUN

#ifdef SINGLE_THREADED
#undef CONCATENATE_INNER
#undef CONCATENATE
#undef ADDSUFFIXTO
#undef INTCONST

#undef INTOF
#undef INTOFCONST
#undef FPOF

#undef SIGN
#undef ABS

#undef COMPUTE_GLOBAL_PARAMS
#undef UPDATE_LOCAL_PARAMS
#undef FPUNION
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

#undef INITRAND_SINGLE
#undef INITRAND_MULTI
#undef GENRAND
#undef GENBIT
#endif /* #ifdef SINGLE_THREADED */

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
