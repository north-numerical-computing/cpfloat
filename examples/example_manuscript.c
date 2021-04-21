/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */
/* SPDX-License-Identifier: LGPL-2.1-or-later                         */

#include <stdio.h>
#include "cpfloat_binary64.h"

int main ()
{
  // Allocate the data structure for target formats and rounding parameters.
  static optstruct *fpopts;
  fpopts = malloc(sizeof(optstruct));

  // Set up the parameters for binary16 target format.
  fpopts->precision = 11;                 // Bits in the significand + 1.
  fpopts->emax = 15;                      // The maximum exponent value.
  fpopts->subnormal = CPFLOAT_SUBN_USE;   // Support for subnormals is on.
  fpopts ->round = CPFLOAT_RND_TP;        // Round toward +infinity.
  fpopts->flip = CPFLOAT_NO_SOFTERR;      // Bit flips are off.
  fpopts->p = 0;                          // Bit flip probability (not used).
  fpopts->explim = CPFLOAT_EXPRANGE_TARG; // Limited exponent in target format.

  // Validate the parameters in fpopts.
  int retval = cpfloat_validate_optstruct(fpopts);
  printf("The validation function returned %d.\n", retval);

  // Initialize a 2x2 matrix with four arbitrary elements
  double X[4] = { (double)1/3, M_PI, M_E, M_SQRT2 };
  double Y[4];
  printf("Values in binary64:\n %.15e %.15e\n %.15e %.15e \n",
         X[0], X[1], X[2], X[3]);

  // Round the values of X to the binary16 format and store in Y
  cpfloat(&Y[0], &X[0], 4, fpopts);
  printf("Rounded to binary16:\n %.15e %.15e\n %.15e %.15e \n",
         Y[0], Y[1], Y[2], Y[3]);

  // Set the precision of the significand to 8 bits,
  // and the maximum exponent to 127, which gives the bfloat16 format
  fpopts ->precision = 8;
  fpopts ->emax = 127;

  // Round the values of X to the bfloat16 and store in Y
  cpfloat(&Y[0], &X[0], 4, fpopts);
  printf("Rounded to bfloat16:\n %.15e %.15e\n %.15e %.15e \n",
         Y[0], Y[1], Y[2], Y[3]);

  return 0;
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
