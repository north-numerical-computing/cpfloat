/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

/*
 * This file is part of CPFLoat.
 *
 * The macros and functions defined here require a number of macros be defined
 * before this header is included. This file should never be included directly
 * in production code. In order to use CPFloat, it suffices to include
 * cpfloat_binary32.h and cpfloat_binary64.h. These two headers include this
 * file multiple times handling the definition of macros correctly.
 */

#ifndef CONCATENATE_INNER
#define CONCATENATE_INNER(arg1, arg2) arg1 ## arg2
#define CONCATENATE(arg1, arg2) CONCATENATE_INNER(arg1, arg2)
#define ADDSUFFIXTO(x) CONCATENATE(x, FUNSUFFIX)

#define MAINFUNNAME cpfloat
#define MASKFUNNAME ADDSUFFIXTO(set_subnormal_masks_)
#define STRUCTNAME ADDSUFFIXTO(fpint)

typedef union {
  INTTYPE intval;
  FPTYPE fpval;
} STRUCTNAME;

#define INTCONST(x) CONCATENATE(x, INTSUFFIX)

#define INTOF(x)(((STRUCTNAME *)(x))->intval)
#define FPOF(x)(((STRUCTNAME)((INTTYPE)(x))).fpval)

#define SIGN(x)(SIGNMASK & INTOF(x))
#define ABS(x)(FPOF((INTTYPE)(ABSMASK & INTOF(x))))

#ifdef PCG_VARIANTS_H_INCLUDED
#define BITTYPE unsigned int
#define BITSEEDTYPE pcg32_random_t
#define INITBIT_SINGLE(seed)                    \
  pcg32_srandom_r(&seed,                        \
                  time(NULL),                   \
                  (intptr_t)&seed)
#define INITBIT_MULTI(seed)                                   \
  pcg32_srandom_r(&seed,                                       \
                  omp_get_thread_num() * 13254 + time(NULL),   \
                  (intptr_t)&seed)
#define GENBIT(seed) (pcg32_random_r(&seed) & (1U << 31))
#else /* #ifdef PCG_VARIANTS_H_INCLUDED */
#define BITTYPE unsigned int
#define BITSEEDTYPE unsigned int
#define INITBIT_SINGLE(seed) seed = time(NULL)
#define INITBIT_MULTI(seed) seed = omp_get_thread_num() * 13254 + time(NULL)
#define GENBIT(seed) (rand_r(&seed) & (1U << 30))
#endif /* #ifdef PCG_VARIANTS_H_INCLUDED */

#define VALFUNNAME CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _validate_optstruct)

static inline
int VALFUNNAME (const optstruct *fpopts) {
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
  size_t i;
  for (i=0; i<nformats; i++) {
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
  if (fpopts->round < -1 || fpopts->round > 9)
    retval = -4;

  // Return 5 if p is required but is not a valid probability.
  if (fpopts->flip && (fpopts->p > 1 || fpopts->p < 0))
    return 5;

  // Return 0 or warning value.
  return retval;

}

static inline
void MASKFUNNAME(const FPTYPE *A,
                 int i,
                 size_t prec,
                 FPTYPE xmin,
                 int emin,
                 INTTYPE leadmask,
                 INTTYPE leadmaskwos,
                 INTTYPE trailmask,
                 size_t *locprec,
                 INTTYPE *locleadmask,
                 INTTYPE *locleadmaskwos,
                 INTTYPE *loctrailmask) {
  // Take care of subnormal values
  if (ABS(A+i) < xmin && emin > 1-DEFEMAX) {
    *locprec = prec - emin - DEFEMAX +
      ((EXPMASK & INTOF(A+i)) >> (DEFPREC - 1));
    *locleadmask = leadmask;
    *locleadmask ^=
      ((INTCONST(1) << (prec - *locprec)) - INTCONST(1)) << (DEFPREC-prec);
    *locleadmaskwos = *locleadmask & ABSMASK;
    *loctrailmask = *locleadmaskwos ^ ABSMASK;
  } else {
    *locprec = prec;
    *locleadmask = leadmask;
    *locleadmaskwos = leadmaskwos;
    *loctrailmask = trailmask;
  }
}
#endif /* #ifndef CONCATENATE_INNER */

#ifdef SINGLE_THREADED
#define FUNNAME CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _sequential)
#else /* #ifdef SINGLE_THREADED */
#define FUNNAME CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _parallel)
#define USE_OPENMP
#endif /* #ifdef SINGLE_THREADED */

/*  Rounding function */
static inline
int FUNNAME(FPTYPE *X,
            const FPTYPE *A,
            const size_t numelem,
            const optstruct *fpopts) {

  int retval = 0;

  size_t prec = fpopts->precision;
  size_t emax = fpopts->explim ? fpopts->emax : DEFEMAX;
  const signed char round = fpopts->round;

  if (prec > DEFPREC) {
    prec = DEFPREC;
    retval = 1;
  }
  if (emax > DEFEMAX) {
    emax = DEFEMAX;
    retval = 2;
  }

  // Derived floating-point parameters.
  const int emin = 1-emax;
  const int esmin = emin-prec+1;
  const FPTYPE xmin = ldexp(1., emin);  // Smallest pos. normal.
  const FPTYPE smin = ldexp(1., esmin); // Smallest pos. subnormal.
  const FPTYPE xmax = ldexp(1., emax) * (2-ldexp(1., 1-prec));

  // A point halfway between the maximum representable number and infinity.
  const FPTYPE xbnd = ldexp(1., emax) * (2-ldexp(1., -prec));

  const FPTYPE ftzthreshold =
    fpopts->subnormal ? smin : xmin; // Flush-to-zero barrier.

  const INTTYPE leadmask = FULLMASK << (DEFPREC-prec); // Bits to keep
  const INTTYPE leadmaskwos = leadmask & ABSMASK; // Bits to keep without sign.
  const INTTYPE trailmask = leadmaskwos ^ ABSMASK; // Bits to discard.

  INTTYPE rndbuf;
  BITTYPE randombit;

  size_t i, locprec;
  INTTYPE locleadmask, locleadmaskwos, loctrailmask;
  #ifdef USE_OPENMP
  #pragma omp parallel                                                  \
    private (locprec, locleadmask, locleadmaskwos, loctrailmask, rndbuf) \
    shared(A, X, fpopts)
  #endif /* #ifdef USE_OPENMP */
  {
    SEEDTYPE seed;
    BITSEEDTYPE bitseed;

    #ifdef USE_OPENMP
    INITRAND_MULTI(seed);
    INITBIT_MULTI(bitseed);
    #else /* #ifdef USE_OPENMP */
    INITRAND_SINGLE(seed);
    INITBIT_SINGLE(bitseed);
    #endif /* #ifdef USE_OPENMP */
    switch (round) {
    case -1: // round-to-zero with ties-to-away
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        if (ABS(A+i) < ftzthreshold) {
          if (ABS(A+i) < ftzthreshold/2)
            X[i] = 0;
          else
            X[i] = FPOF(INTOF(&ftzthreshold) | SIGN(A+i));
        } else {
          X[i] = FPOF(((INTOF(A+i) & ABSMASK)
                       + (INTCONST(1) << (DEFPREC-1-locprec))) & locleadmask);
          // Overflow.
          if (X[i] >= xbnd)
            X[i] = INFINITY;
          // Restore sign.
          X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
        }
      }
      break;

    case 0: // round-to-zero with ties-to-zero
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        if (ABS(A+i) < ftzthreshold) {
          if (ABS(A+i) <= ftzthreshold/2)
            X[i] = 0;
          else
            X[i] = FPOF(INTOF(&ftzthreshold) | SIGN(A+i));
        } else {
          X[i] = FPOF(((INTOF(A+i) & ABSMASK)
                       + (loctrailmask>>1)) & locleadmask);
          // Overflow.
          if (ABS(A+i) >= xbnd)
            X[i] = INFINITY;
          // Restore sign.
          X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
        }
      }
      break;

    case 1: // round-to-zero with ties-to-even
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        if (ABS(A+i) < ftzthreshold) {
          if (ABS(A+i) < ftzthreshold/2
              || (ABS(A+i) == ftzthreshold/2 && fpopts->subnormal))
            X[i] = 0;
          else
            X[i] = FPOF(INTOF(&ftzthreshold));
        } else if (ABS(A+i) >= xbnd) {
          X[i] = INFINITY;
        } else {
          INTTYPE absval = INTOF(A+i) & ABSMASK;
          INTTYPE LSB = ((absval >> (DEFPREC-locprec)) & INTCONST(1))
            | (locprec == 1 && DEFEMAX != emax); // Hidden bit is one.
          X[i] = FPOF((absval + ((loctrailmask >> 1) + LSB)) & locleadmask);
        }
        // Restore sign.
        X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
      }
      break;

    case 2: // round-toward-positive
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        if (ABS(A+i) < ftzthreshold) {
          X[i] = A[i] > 0 ? ftzthreshold : 0;
        } else {
          X[i] = FPOF(INTOF(A+i) & locleadmask);
          if (SIGN(A+i) == 0) // Add ulp if x is positive.
            X[i] = FPOF((INTOF(A+i) + loctrailmask) & locleadmask);
          // Overflow.
          if (X[i] > xmax)
            X[i] = INFINITY;
          else if (X[i] < -xmax && X[i] != - INFINITY)
            X[i] = -xmax;
        }
      }
      break;

    case 3: // round-toward-negative
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        if (ABS(A+i) < ftzthreshold)
          X[i] = A[i] >= 0 ? 0 : -ftzthreshold;
        else {
          X[i] = FPOF(INTOF(A+i) & locleadmask);
          if (SIGN(A+i)) // Subtract ulp if x is positive.
            X[i] = FPOF((INTOF(A+i) + loctrailmask) & locleadmask);
          // Overflow.
          if (X[i] > xmax && X[i] != INFINITY)
            X[i] = xmax;
          else if (X[i] < -xmax)
            X[i] = -INFINITY;
        }
      }
      break;

    case 4: // round-toward-zero
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        X[i] = FPOF(INTOF(A+i) & locleadmaskwos);
        // Underflow and overflow.
        if (ABS(A+i) < ftzthreshold)
          X[i] = 0;
        else {
          if (X[i] > xmax && X[i] != INFINITY)
            X[i] = xmax;
          // Restore sign.
          X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
        }
      }
      break;

    case 5: // stochastic rounding with proportional probabilities
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        // Underflow and overflow.
        X[i] = ABS(A+i);
        if (X[i] < ftzthreshold){
          int expx = ((EXPMASK & INTOF(A+i)) >> (DEFPREC - 1)) - DEFEMAX;
          INTTYPE trailfrac = (INTOF(A+i) & (FULLMASK >> NLEADBITS));
          if (expx == -DEFEMAX) { // No implicit bit (x is subnormal).
            expx += 1;
          } else { // Make implicit bit explicit (x is normal).
            trailfrac = (INTOF(A+i) & (FULLMASK >> NLEADBITS))
              | (INTCONST(1) << (DEFPREC-1));
          }
          int localemin = fpopts->subnormal ? emin - (int)prec + 1: emin;
          int expdiff = localemin - expx;
          INTTYPE rnd = GENRAND(seed) >> 1;
          // Shift fraction of A[i] left or right as needed.
          if (expdiff <= NLEADBITS - 1)
            trailfrac <<= NLEADBITS - 1 - expdiff;
          else
            trailfrac >>= expdiff - (NLEADBITS - 1);
          if (trailfrac > rnd) {
            X[i] = ftzthreshold;
          } else {
            X[i] = 0;
            continue;
          }
        } else if (X[i] < xbnd) {
          rndbuf = GENRAND(seed) & loctrailmask;
          if (X[i] > xmax) {
            rndbuf = rndbuf >> 1;
            locleadmask = locleadmask >> 1;
          }
          X[i] = FPOF(((INTOF(A+i) & ABSMASK) + rndbuf) & locleadmask);
        }
        //Overflow.
        if (X[i] >= xbnd)
          X[i] = INFINITY;
        // Restore sign.
        X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
      }
      break;

    case 6: // stochastic rounding with equal probabilities
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        // Underflow and overflow.
        if (ABS(A+i) < ftzthreshold && A[i] != 0) {
          randombit = GENBIT(bitseed);
          X[i] = randombit ? ftzthreshold : 0;
          X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
        } else if (ABS(A+i) > xmax && ABS(A+i) != INFINITY) {
          randombit = GENBIT(bitseed);
          X[i] = randombit ? INFINITY  : xmax;
          X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
        } else if ((INTOF(A+i) & loctrailmask)) { // Not exactly representable.
          randombit = GENBIT(bitseed);
          X[i] = FPOF(INTOF(A+i) & locleadmaskwos);
          if (randombit)
            X[i] = FPOF(INTOF(X+i) + (INTCONST(1) << (DEFPREC-locprec)));
          // Restore sign.
          X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
        } else
          X[i] = A[i];
      }
      break;

    case 7: // round-to-odd
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        MASKFUNNAME(A, i, prec, xmin, emin,
                    leadmask, leadmaskwos, trailmask,
                    &locprec, &locleadmask,
                    &locleadmaskwos, &loctrailmask);
        // Underflow.
        if (ABS(A+i) < ftzthreshold && A[i] != 0)
          X[i] = ftzthreshold;
        else {
          X[i] = FPOF(INTOF(A+i) & locleadmaskwos);
          if ((loctrailmask & INTOF(A+i)) // Not exactly representable.
              && (locprec > 1)) // Set the last bit to one if stored explicitly.
            X[i] = FPOF(INTOF(X+i) | (INTCONST(1) << (DEFPREC-locprec)));
          // Overflow.
          if (X[i] > xmax && X[i] != INFINITY)
            X[i] = xmax;
        }
        // Restore sign.
        X[i] = FPOF(INTOF(X+i) | SIGN(A+i));
      }
      break;

    default: // No rounding if unknown mode specified.
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        X[i] = A[i];
      }
      break;
    }

    // Bit flips.
    if (fpopts->flip) {
      #ifdef USE_OPENMP
      #pragma omp for
      #endif /* #ifdef USE_OPENMP */
      for (i=0; i<numelem; i++){
        if (rand() / (FPTYPE)RAND_MAX < fpopts->p) {
          X[i] = FPOF(INTOF(X+i) ^ (INTCONST(1) << rand() % (DEFPREC - 1)));
        }
      }
    }
  }

  return retval;
}

#ifdef _OPENMP
#ifndef USE_OPENMP
#define FUNNAME_SINGLE CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _sequential)
#define FUNNAME_MULTI  CONCATENATE(ADDSUFFIXTO(MAINFUNNAME), _parallel)
#define FUNNAME_COMBO  ADDSUFFIXTO(MAINFUNNAME)
int FUNNAME_COMBO(FPTYPE *X,
                  const FPTYPE *A,
                  const size_t numelem,
                  const optstruct *fpopts) {
  if (numelem < CONCATENATE(OPENMP_THRESHOLD_, FPTYPE))
    return FUNNAME_SINGLE(X, A, numelem, fpopts);
  else
    return FUNNAME_MULTI(X, A, numelem, fpopts);
}
#undef FUNNAME_SINGLE
#undef FUNNAME_MULTI
#undef FUNNAME_COMBO
#endif /*#ifdef USE_OPENMP */
#else /* #ifdef _OPENMP */
int ADDSUFFIXTO(MAINFUNNAME)(FPTYPE *X,
                             const FPTYPE *A,
                             const size_t numelem,
                             const optstruct *fpopts) {
  return FUNNAME(X, A, numelem, fpopts);
}
#endif /* #ifdef _OPENMP */

#undef FUNNAME

#ifdef USE_OPENMP
#undef USE_OPENMP
#endif /* #ifdef USE_OPENMP */

#ifdef SINGLE_THREADED
#undef CONCATENATE_INNER
#undef CONCATENATE
#undef ADDSUFFIXTO
#undef INTCONST

#undef INTOF
#undef FPOP

#undef SIGN
#undef ABS

#undef MASKFUNNAME
#undef STRUCTNAME
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

#undef SEEDTYPE
#undef INITRAND_SINGLE
#undef INITRAND_MULTI
#undef GENRAND
#undef GENBIT
#undef NRNDBITS
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
