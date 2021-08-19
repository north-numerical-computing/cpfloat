% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

function retval = cpfloat_compile(varargin)
%COMPILE_CHOPFAST    Compile MEX interface to the CPFloat Library.
%   COMPILE_CHOPFAST() compiles the MEX function chopfast using the default C
%   compiler. The function expects all the header files of the CPFloat Library
%   as well as the file pcg_variants.h from the PCG Library to be in the current
%   working directory. The function attempts to use the OpenMP library, if
%   available.
%
%   COMPILE_CHOPFAST('cpfloatdir',CPFLOATDIR) looks for the header files of the
%   CPFloat Library in CPFLOATDIR rather than in the current working directory.
%
%   COMPILE_CHOPFAST('pcgpath',PCGPATH) sets the root directory of the PCG
%   random number generator to PCGPATH instead of ./pcg-c/.
%
%   COMPILE_CHOPFAST('pcgvariants',PCGVARIANTS) specifies that the path of the
%   header file pcg_variants.h is PCGVARIANTS. The defauls value is
%   PCGPATH/include/pcg_variants.h.
%
%   COMPILE_CHOPFAST('pcglib',PCGLIB) specifies that the path the static library
%   libpcg_random.a is PCGLIB. The defauls value is PCGPATH/src/libpcg_random.a.
%
%   COMPILE_CHOPFAST('compilerpath',COMPILERPATH) uses the compiler COMPILERPATH
%   instead of the default C compiler.

  retval = true;

  p = inputParser;
  addParameter(p, 'cpfloatdir', '', @ischar);
  addParameter(p, 'pcgpath', './pcg-c/', @ischar);
  addParameter(p, 'pcgvariants', '', @ischar);
  addParameter(p, 'pcglib', '', @ischar);
  addParameter(p, 'compilerpath', '', @ischar);
  parse(p,varargin{:});
  cpfloatdir = p.Results.cpfloatdir;
  pcgpath = p.Results.pcgpath;
  pcgvariants = p.Results.pcgvariants;
  pcglib = p.Results.pcglib;
  compilerpath = p.Results.compilerpath;

  coptions = '-std=gnu99 -O3 -march=native';

  % Try to find the PCG library.
  if (isempty(pcgvariants))
    pcgvariants = sprintf('%s/include/pcg_variants.h', pcgpath);
  end
  if (isempty(pcglib))
    pcglib = sprintf('%s/src/libpcg_random.a', pcgpath);
  end
  if exist(pcgvariants, 'file') && exist(pcglib, 'file')
    coptions = sprintf('%s -include%s', coptions, pcgvariants);
    clibs = sprintf('-L%s/', fileparts(pcglib));
  else
    pcglib = '';
    clibs = '';
  end
  if ~isempty(cpfloatdir)
    coptions = sprintf('-I%s %s', cpfloatdir, coptions);
  end

  usingoctave = exist('OCTAVE_VERSION', 'builtin');
  if usingoctave
    if ~isempty(compilerpath)
      setenv("CC", compilerpath);
      setenv("CXX", compilerpath);
      setenv("DL_LD", compilerpath);
    end
    setenv("CFLAGS", sprintf("-fopenmp %s", coptions));
    setenv("LDFLAGS", sprintf("-fopenmp %s", clibs));
    if isempty(pcglib)
      [output, status] = mkoctfile('cpfloat.c', '--mex', '--verbose');
    else
      [output, status] = mkoctfile('cpfloat.c', pcglib, '--mex', '--verbose');
    end
    if status ~= 0
      warning('Compilation error, trying to compile without OpenMP.');
      retval = false;
      setenv("CFLAGS", coptions);
      setenv("LDFLAGS", clibs)
      if isempty(pcglib)
        [output, status] = mkoctfile('cpfloat.c', '--mex', '--verbose');
      else
        [output, status] = mkoctfile('cpfloat.c', pcglib, '--mex', '--verbose');
      end
    end
  else
    if isempty(compilerpath)
      compiler_string = '';
    else
      compiler_string = ['CC="' compilerpath '"'];
    end
    try
      mex('cpfloat.c', pcglib, '-silent',...
          compiler_string,...
          [sprintf('CFLAGS=$CFLAGS %s -fopenmp ', coptions)],...
          [sprintf('LDFLAGS=$LDFLAGS %s -fopenmp ', clibs)]);
    catch
      warning('Compilation error, trying to compile without OpenMP.');
      retval = false;
      mex('cpfloat.c', pcglib, '-silent',...
          compiler_string,...
          [sprintf('CFLAGS="$CFLAGS %s" ', coptions)],...
          [sprintf('LDFLAGS="$LDFLAGS %s" ', clibs)]);
    end
  end
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
