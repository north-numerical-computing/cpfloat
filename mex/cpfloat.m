%CPFLOAT    Round floating point numbers to lower precision.
%
%   [Y, OPTIONS] = CPFLOAT(X, FPOPTS) returns a matrix Y containing the elements
%   of X rounded to a lower-precision floating-point format (the target format).
%   The function can be used to simulate the occurrence of soft errors in the
%   rounded values. X must be a real matrix with entries of class 'single' or
%   'double' (the storage format), and the output matrix Y will be a real matrix
%   of the same size with entries of the same class. The parameters that
%   describe the target format, the rounding mode, and the likelihood of soft
%   errors are stored by the function in persistent memory, and are preserved
%   across multiple calls to CPFLOAT. The internal configuration can be modified
%   by means of the structure FPOPTS, whose fields are discussed in detail
%   below. The parameters for which a new configuration value is not specified
%   take the default value on the first invocation of CPFLOAT, and keep their
%   previous values on subsequent calls. The parameters of the current
%   configuration are returned in the second output argument OPTIONS, a
%   structure with the same fields as FPOPTS.
%
%   The fields of FPOPTS are interpreted as follows.
%
%   * The string FPOPTS.format specifies the target floating-point format.
%     Possible values are:
%       'b', 'bf16', 'bfloat16'            for Intel bfloat16;
%       'h', 'fp16', 'binary16', 'half'    for IEEE binary16 (half precision);
%       't', 'tf32', 'TeensorFloat-32'     for NVIDIA TensorFloat-32;
%       's', 'fp32', 'binary32', 'single'  for IEEE binary32 (single precision);
%       'd', 'fp64', 'binary64', 'double'  for IEEE binary64 (double precision);
%       'c', 'custom'                      for a custom-precision format.
%     In order to use a custom precision, the parameters of the floating-point
%     format must be supplied using the FPOPTS.params field. The default value
%     for this field is 'h'.
%
%   * The two-element vector FPOPTS.params specifies the parameters of the
%     target floating-point format, and is ignored unless FPOPTS.format is set
%     to either 'c' or 'custom'. The vector has the form [PRECISION , EMAX],
%     where PRECISION and EMAX are positive integer representing the number of
%     binary digits in the fraction and the maximum exponent of the target
%     format, respectively. The minimum exponent is assumed to be 1 - EMAX. The
%     default value of this field is the vector [11, 15].
%
%   * The scalar FPOPTS.subnormal specifies the support for subnormal numbers.
%     The target floating-point format will not support subnormal numbers if
%     this field is set to 0, and will support them otherwise. The default value
%     for this field is 0 if the target format is 'bfloat16' and 1 otherwise.
%
%   * The scalar FPOPTS.explim specifies the support for an extended exponent
%     range. The target floating-point format will have the exponent range of
%     the storage format ('single' or 'double', depending on the class of X) if
%     this field is set to 0, and the exponent range of the format specified in
%     FPOPTS.format otherwise. The default value for this field is 1.
%
%   * The scalar FPOPTS.round specifies the rounding mode. Possible values are:
%       -1 for round-to-nearest with ties-to-away;
%        0 for round-to-nearest with ties-to-zero;
%        1 for round-to-nearest with ties-to-even;
%        2 for round-toward-plus-infinity;
%        3 for round-toward-minus-infinity;
%        4 for round-toward-zero;
%        5 for round-stochastic with proportional probabilities;
%        6 for round-stochastic with equal probabilities;
%        7 for round-to-odd;
%        8 for no rounding.
%      Any other value results in no rounding. The default value for this field
%      is 1.
%
%   * The scalar FPOPTS.flip specifies whether the function should simulate the
%     occurrence of soft errors striking the elements of Y. If this field is not
%     set to 0, then the fraction of each element of Y will have a randomly
%     chosen bit flipped with probability FPOPTS.p. If the exponent range of the
%     storage format is larger than that of the target format, then subnormal
%     numbers might be stored as normal numbers, in which case the bit flip
%     cannot strike the leading bit of the representation. The default value for
%     this field is 0.
%
%   * The scalar FPOPTS.p specifies the probability of bit flips. If FPOPTS.flip
%     is not set to zero, then the value of this field must be a valid
%     probability, that is, a real number in the interval [0,1]. The default
%     value for this field is 0.5.
%
%   The interface of CPFLOAT is fully compatible with that of the MATLAB
%   function CHOP available at https://github.com/higham/chop.

% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

% CPFloat - Custom Precision Floating-point numbers.
%
% Copyright 2020 Massimiliano Fasi and Mantas Mikaitis
%
% This library is free software; you can redistribute it and/or modify it under
% the terms of the GNU Lesser General Public License as published by the Free
% Software Foundation; either version 2.1 of the License, or (at your option)
% any later version.
%
% This library is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
% details.
%
% You should have received a copy of the GNU Lesser General Public License along
% with this library; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
