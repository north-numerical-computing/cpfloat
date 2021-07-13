/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "mex.h"
#include "cpfloat_binary32.h"
#include "cpfloat_binary64.h"

static optstruct *fpopts;
void clearfpopts() {
  if (fpopts != NULL)
    mxFree(fpopts);
}

/********************
 * GATEWAY FUNCTION *
 ********************/
void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[]) {

  /* Check for correct number of arguments. */
  if(nrhs > 3) {
    mexErrMsgIdAndTxt("cpfloat:nrhs",
                      "Chopfast requires at most three input arguments");
  }

  /* Allocate fpopts struct and set fields to default. */
  if (fpopts == NULL) {
    fpopts = mxCalloc(1, sizeof(optstruct));
    mexMakeMemoryPersistent(fpopts);
    mexAtExit(clearfpopts);

    strcpy(fpopts->format, "h");
    fpopts->precision = 11; // t
    fpopts->emax = 15; // emax
    fpopts->subnormal = CPFLOAT_SUBN_USE;
    fpopts->explim = CPFLOAT_EXPRANGE_TARG;
    fpopts->round = CPFLOAT_RND_NE;
    fpopts->flip = CPFLOAT_NO_SOFTERR;
    fpopts->p = 0.5;
    fpopts->bitseed = NULL;
    fpopts->randseedf = NULL;
    fpopts->randseed = NULL;
  }

  /* Parse second argument and populate fpopts structure. */
  if (nrhs > 1) {
    if(!mxIsStruct(prhs[1])) {
      mexErrMsgIdAndTxt("cpfloat:invalidstruct",
                        "Second argument must be a struct.");
    } else {
      mxArray *tmp = mxGetField(prhs[1], 0, "format");

      if (tmp != NULL) {
        if (mxGetM(tmp) == 0 && mxGetN(tmp) == 0)
          /* Use default format, for compatibility with chop. */
          strcpy(fpopts->format, "h");
        else if (mxGetClassID(tmp) == mxCHAR_CLASS)
          strcpy(fpopts->format, mxArrayToString(tmp));
      }
      tmp = mxGetField(prhs[1], 0, "params");
      if ((tmp != NULL) &&
          (strcmp(fpopts->format, "c")
           && strcmp(fpopts->format, "custom")))
        mexWarnMsgIdAndTxt("cpfloat:ignoredparams",
                           "Floating-point parameters ignored.");
      /* Populate fpopts->params according to fpopts->format. */
      if (!strcmp(fpopts->format, "b") ||
          !strcmp(fpopts->format, "bfloat16") ||
          !strcmp(fpopts->format, "bf16")) {
        fpopts->precision = 8;   // t
        fpopts->emax = 127;      // emax
        fpopts->subnormal = CPFLOAT_SUBN_RND; // default for bfloat16
      } else if (!strcmp(fpopts->format, "h") ||
                 !strcmp(fpopts->format, "half") ||
                 !strcmp(fpopts->format, "binary16") ||
                 !strcmp(fpopts->format, "fp16")) {
        fpopts->precision = 11; // t
        fpopts->emax = 15;      // emax
      } else if (!strcmp(fpopts->format, "t") ||
                 !strcmp(fpopts->format, "TensorFloat-32") ||
                 !strcmp(fpopts->format, "tf32")) {
        fpopts->precision = 11;  // t
        fpopts->emax = 127;      // emax
      } else if (!strcmp(fpopts->format, "s") ||
                 !strcmp(fpopts->format, "single") ||
                 !strcmp(fpopts->format, "binary32") ||
                 !strcmp(fpopts->format, "fp32")) {
        fpopts->precision =  24; // t
        fpopts->emax = 127;      // emax
      } else if (!strcmp(fpopts->format, "d") ||
                 !strcmp(fpopts->format, "double") ||
                 !strcmp(fpopts->format, "binary64") ||
                 !strcmp(fpopts->format, "fp64")) {
        fpopts->precision =   53; // t
        fpopts->emax = 1023;      // emax
      } else if (!strcmp(fpopts->format, "c") ||
                 !strcmp(fpopts->format, "custom")) {
        if ((tmp != NULL) && (mxGetClassID(tmp) == mxDOUBLE_CLASS)) {
          fpopts->precision = ((double *)mxGetData(tmp))[0]; // t
          fpopts->emax = ((double *)mxGetData(tmp))[1];      // emax
        } else {
          mexErrMsgIdAndTxt("cpfloat:invalidparams",
                            "Invalid floating-point parameters specified.");
        }
      } else {
        mexErrMsgIdAndTxt("cpfloat:invalidformat",
                          "Invalid floating-point format specified.");
      }
      /* Set default values o be compatible with MATLAB chop. */
      tmp = mxGetField(prhs[1], 0, "subnormal");
      if (tmp != NULL) {
        if (mxGetM(tmp) == 0 && mxGetN(tmp) == 0)
          fpopts->subnormal = CPFLOAT_SUBN_USE;
        else if (mxGetClassID(tmp) == mxDOUBLE_CLASS)
          fpopts->subnormal = *((double *)mxGetData(tmp));
      }
      tmp = mxGetField(prhs[1], 0, "explim");
      if (tmp != NULL) {
        if (mxGetM(tmp) == 0 && mxGetN(tmp) == 0)
          fpopts->explim = 1;
        else if (mxGetClassID(tmp) == mxDOUBLE_CLASS)
          fpopts->explim = *((double *)mxGetData(tmp));
      }
      tmp = mxGetField(prhs[1], 0, "round");
      if (tmp != NULL) {
        if (mxGetM(tmp) == 0 && mxGetN(tmp) == 0)
          fpopts->round = CPFLOAT_RND_NE;
        else if (mxGetClassID(tmp) == mxDOUBLE_CLASS)
          fpopts->round = *((double *)mxGetData(tmp));
      }
      tmp = mxGetField(prhs[1], 0, "flip");
      if (tmp != NULL) {
        if (mxGetM(tmp) == 0 && mxGetN(tmp) == 0)
          fpopts->flip = CPFLOAT_NO_SOFTERR;
        else if (mxGetClassID(tmp) == mxDOUBLE_CLASS)
          fpopts->flip = *((double *)mxGetData(tmp));
      }
      tmp = mxGetField(prhs[1], 0, "p");
      if (tmp != NULL) {
        if (mxGetM(tmp) == 0 && mxGetN(tmp) == 0)
          fpopts->p = 0.5;
        else if (mxGetClassID(tmp) == mxDOUBLE_CLASS)
          fpopts->p = *((double *)mxGetData(tmp));
      }
    }
  }

  int algorithm;
  if (nrhs > 2) {
    double *tmp = (double *)mxGetData(prhs[2]);
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
        || *tmp != round(*tmp)
        || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 1)
      mexErrMsgIdAndTxt("cpfloat:invalidalgorithm",
                        "Third parameters must be an integer.");
    algorithm = (int)(*tmp);
    #ifdef _OPENMP
    if (algorithm < -omp_get_max_threads())
      algorithm = -omp_get_max_threads();
    else if (algorithm > omp_get_max_threads())
      algorithm = omp_get_max_threads();
    #endif
  } else
    algorithm = 0;

  /* Parse first argument. */
  if (nrhs > 0) {
    if (!mxIsNumeric(prhs[0])
        || (!mxIsDouble(prhs[0]) && !mxIsSingle(prhs[0]))
        || mxIsComplex(prhs[0])
        || (mxGetNumberOfDimensions(prhs[0]) != 2)) {
      mexErrMsgIdAndTxt("cpfloat:invalidmatrix",
                        "First argument must be a 2D real numeric array.");
    }

    mwSize maxfbits, maxebits;
    if (mxIsSingle(prhs[0])) {
      if (!strcmp(fpopts->format, "d") ||
          !strcmp(fpopts->format, "double") ||
          !strcmp(fpopts->format, "binary64") ||
          !strcmp(fpopts->format, "fp64")) {
        mexErrMsgIdAndTxt("cpfloat:invalidformat",
                          "Target format is too large.");
      } else {
        maxfbits = fpopts->round<=1 ? 11 : 23;
        maxebits = 127;
      }
    } else if(mxIsDouble(prhs[0])) {
      maxfbits = fpopts->round<=1 ? 25 : 52;
      maxebits = 1023;
    }
    if (fpopts->precision > maxfbits || fpopts->emax > maxebits)
      if (!strcmp(fpopts->format, "c") || !strcmp(fpopts->format, "custom"))
        mexErrMsgIdAndTxt("cpfloat:invalidparams",
                          "Invalid floating-point parameters selected.");

    /* Allocate and compute first output. */
    mwSize m, n;
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    mwSize dims[2];
    dims[0] = m;
    dims[1] = n;

    if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
      double *A = (double *)mxGetData(prhs[0]);
      plhs[0] = mxCreateNumericArray(2, dims,mxDOUBLE_CLASS, mxREAL);
      double *X = (double *)mxGetData(plhs[0]);
      #ifdef _OPENMP
      if (algorithm == 0) {
        cpfloat(X, A, m*n, fpopts);
      } else if (algorithm == 1){
        cpfloat_sequential(X, A, m*n, fpopts);
      } else if (algorithm > 0) {
        omp_set_num_threads(algorithm);
        cpfloat(X, A, m*n, fpopts);
      } else {
        omp_set_num_threads(-algorithm);
        cpfloat_parallel(X, A, m*n, fpopts);
      }
      #else
      cpfloat(X, A, m*n, fpopts);
      #endif
    } else if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS) {
      float *A = (float *)mxGetData(prhs[0]);
      plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS,mxREAL);
      float *X = (float *)mxGetData(plhs[0]);
      #ifdef _OPENMP
      if (algorithm == 0) {
        cpfloatf(X, A, m*n, fpopts);
      } else if (algorithm == 1){
        cpfloatf_sequential(X, A, m*n, fpopts);
      } else if (algorithm > 0) {
        omp_set_num_threads(algorithm);
        cpfloatf(X, A, m*n, fpopts);
      } else {
        omp_set_num_threads(-algorithm);
        cpfloatf_parallel(X, A, m*n, fpopts);
      }
      #else
      cpfloatf(X, A, m*n, fpopts);
      #endif
    } else {
      mexErrMsgIdAndTxt("cpfloat:invalidmatrix",
                        "First argument must be a numeric array.");
    }
  } else {
    mwSize dims[2];
    dims[0] = 0;
    dims[1] = 0;
    plhs[0] = mxCreateNumericArray(2, dims,mxDOUBLE_CLASS, mxREAL);
  }

  /* Allocate and return second output. */
  if (nlhs > 1) {
    const char* field_names[] = {"format", "params", "subnormal", "round",
                                 "flip", "p", "explim"};
    mwSize dims[2] = {1, 1};
    plhs[1] = mxCreateStructArray(2, dims, 7, field_names);
    mxSetFieldByNumber(plhs[1], 0, 0, mxCreateString(fpopts->format));

    mxArray *outparams = mxCreateDoubleMatrix(1,2,mxREAL);
    double *outparamsptr = mxGetData(outparams);
    outparamsptr[0] = fpopts->precision;
    outparamsptr[1] = fpopts->emax;
    mxSetFieldByNumber(plhs[1], 0, 1, outparams);

    mxArray *outsubnormal = mxCreateDoubleMatrix(1,1,mxREAL);
    double *outsubnormalptr = mxGetData(outsubnormal);
    outsubnormalptr[0] = fpopts->subnormal;
    mxSetFieldByNumber(plhs[1], 0, 2, outsubnormal);

    mxArray *outround = mxCreateDoubleMatrix(1,1,mxREAL);
    double *outroundptr = mxGetData(outround);
    outroundptr[0] = fpopts->round;
    mxSetFieldByNumber(plhs[1], 0, 3, outround);

    mxArray *outflip = mxCreateDoubleMatrix(1,1,mxREAL);
    double *outflipptr = mxGetData(outflip);
    outflipptr[0] = fpopts->flip;
    mxSetFieldByNumber(plhs[1], 0, 4, outflip);

    mxArray *outp = mxCreateDoubleMatrix(1,1,mxREAL);
    double *outpptr = mxGetData(outp);
    outpptr[0] = fpopts->p;
    mxSetFieldByNumber(plhs[1], 0, 5, outp);

    mxArray *outexplim = mxCreateDoubleMatrix(1,1,mxREAL);
    double *outexplimptr = mxGetData(outexplim);
    outexplimptr[0] = fpopts->explim;
    mxSetFieldByNumber(plhs[1], 0, 6, outexplim);

  }
  if (nlhs > 2)
    mexErrMsgIdAndTxt("cpfloat:invalidnargout",
                      "This function returns at most two valaues.");

}

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
