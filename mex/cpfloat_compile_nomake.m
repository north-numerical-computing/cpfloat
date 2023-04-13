% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

% This MATLAB/Octave script attempts to build the MEX interface to CPFloat on
% systems where the make tool is not available. The code provides only minimal
% functionalities, but should produce a MEX file on a machines where the C
% building environment is configured correctly.

% Absolute path of the C compiler to be used to build the MEX interface.
% If the string is left empty, the default C compiler will be used.
compilerpath = '';

% Absolute path of the source code of cpfloat. By default, the script
% assumes that it is being run from the cpfloat/mex/ folder.
cpfloat_dir = fileparts(pwd);

% Compile MEX interface.
cpfloat_srcdir = fullfile(cpfloat_dir, 'src');
retval = cpfloat_compile('cpfloatdir', cpfloat_srcdir,...
                         'compilerpath', compilerpath);

% If parallel compilation was successful, auto-tune the threshold.
if retval
  cpfloat_autotune('cpfloatdir', cpfloat_srcdir);
  cpfloat_compile('cpfloatdir', cpfloat_srcdir,...
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
