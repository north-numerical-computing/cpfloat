% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

% This MATLAB/Octave script attempts to build the MEX interface to CPFloat on
% systems where the make tool is not available. The code provides only minimal
% functionalities, but should produce a MEX file on a machines where the C
% building environment is configured correctly.

% Absolute path of the C compiler to be used to build the MEX interface.
% If the string is left empty, the default C compiler will be used.
compilerpath = '';

% Download and compile the PCG Library.
try
  unzip('https://codeload.github.com/imneme/pcg-c/zip/refs/heads/master');
  cd('pcg-c-master/src/');
  system('make libpcg_random.a');
  cd('../../');
catch
  if ~exist('../include/pcg_variants.h', 'file')
    warning('Unable to downlaod the PCG Library header file.');
    warning('The standard C library RNG will be used.');
  end
end

% Compile MEX interface.
retval = cpfloat_compile('cpfloatdir', '../src/',...
                         'pcgpath', './pcg-c-master/',...
                         'compilerpath', compilerpath);

% If parallel compilation was successful, auto-tune the threshold.
if retval
  cpfloat_autotune('cpfloatdir', '../src/');
  cpfloat_compile('cpfloatdir', '../src/',...
                  'pcgpath', './pcg-c-master/',...
                  'compilerpath', compilerpath);
end

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
