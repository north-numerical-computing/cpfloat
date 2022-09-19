% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

function cpfloat_autotune(varargin)
%CPFLOAT_AUTOTUNE    Autotune MEX interface to the CPFloat Library.
%   CPFLOAT_AUTOTUNE() runs the function CPFLOAT with inputs of class 'single'
%   and 'double' and computes the size at which switching from the sequential to
%   the parallel implementation becomes beneficial. The functions generates the
%   two files cpfloat_threshold_binary32.h and cpfloat_threshold_binary64.h in
%   the current working directory.
%
%   CPFLOAT_AUTOTUNE('cpfloatdir',CPFLOATDIR) places the output files in the
%   folder CPFLOATDIR instead of the current woking directory. CPFLOATDIR must
%   be an existing folder.

  fpopts.format = 'h';
  fpopts.subnormal = 1;
  fpopts.round = 1;
  fpopts.flip = 0;
  fpopts.p = 0.5;
  fpopts.explim = 1;

  p = inputParser;
  addParameter(p, 'cpfloatdir', './', @ischar);
  if exist('maxNumCompThreads', 'builtin')
    addParameter(p, 'nthreads', maxNumCompThreads(), ...
                 @(x)(isscalar(x) && round(x) == x));
  else
    pkg load parallel
    addParameter(p, 'nthreads', parcellfun_set_nproc(Inf), ...
                 @(x)(isscalar(x) && round(x) == x));
  end
  parse(p,varargin{:});
  cpfloatdir = p.Results.cpfloatdir;
  nthreads = p.Results.nthreads;

  ntests = 100;

  fprintf('Test using %d OpenMP threads.\n', nthreads);
  if exist('timeit', 'builtin')
    parfaster = @(n, fpopts, ntests, fpclass)...
        parfaster_timeit(n, fpopts, ntests, nthreads, fpclass);
  else
    parfaster = @(n, fpopts, ntests, fpclass)...
        parfaster_tictoc(n, fpopts, ntests, nthreads, fpclass);
  end

  docstring =[
      '/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */\n',...
      '/* SPDX-License-Identifier: LGPL-2.1-or-later                         */\n',...
      '\n',...
      '/**\n',...
      ' * @file %s_threshold_%s.h\n',...
      ' * @brief Size of smallest `%s` array on which to use',...
      ' multiple OpenMP threads.\n',...
      ' */\n',...
      '\n',...
      '/**\n',...
      ' * @brief Size of smallest array on which %s() uses multiple threads.\n',...
      ' *\n',...
      ' * @details Threshold for switching between %s_sequential() and\n',...
      ' * %s_parallel() in %s(). The value of this constant is ignored\n',...
      ' * if the file that includes cpfloat_%s.h is compiled without OpenMP\n',...
      ' * support.\n',...
      ' */\n'];

  % Binary32
  fpclass = 'single';
  nmin = 1;
  nmax = 1;
  while(~parfaster(nmax, fpopts, ntests, fpclass))
    nmax = nmax * 2;
  end
  nmid = round((nmax + nmin) / 2);
  while(nmid ~= nmin && nmid ~= nmax)
    if(parfaster(nmid, fpopts, ntests, fpclass))
      nmax = nmid;
    else
      nmin = nmid;
    end
    nmid = round((nmax + nmin) / 2);
  end
  filename = sprintf('%s/cpfloat_threshold_binary32.h', cpfloatdir);
  fid = fopen(filename, 'w');
  fprintf(fid, docstring, 'cpfloat', 'binary32', 'float',...
          'cpfloatf', 'cpfloatf', 'cpfloatf', 'cpfloatf', 'binary32');
  fprintf(fid, "#define OPENMP_THRESHOLD_float %d", nmax);
  fclose(fid);

  % Binary64
  nmin = 1;
  nmax = 1;
  while(~parfaster(nmax, fpopts, ntests, 'double'))
    nmax = nmax * 2;
  end
  nmid = round((nmax + nmin) / 2);
  while(nmid ~= nmin && nmid ~= nmax)
    if(parfaster(nmid, fpopts, ntests, fpclass))
      nmax = nmid;
    else
      nmin = nmid;
    end
    nmid = round((nmax + nmin) / 2);
  end
  filename = sprintf('%s/cpfloat_threshold_binary64.h', cpfloatdir);
  fid = fopen(filename, 'w');
  fprintf(fid, docstring, 'cpfloat', 'binary64', 'double',...
          'cpfloat', 'cpfloat', 'cpfloat', 'cpfloat', 'binary64');
  fprintf(fid, "#define OPENMP_THRESHOLD_double %d", nmax);
  fclose(fid);

  function res = parfaster_timeit(n, fpopts, ~, nthreads, fpclass)
    X = rand(n, 1, fpclass);
    funseq = @()(cpfloat(X, fpopts, 1));
    seqtime = timeit(funseq);
    funseq = @()(cpfloat(X, fpopts, -nthreads));
    partime = timeit(funseq);
    res = partime < seqtime;
    fprintf('[%7d]   %.5e   %.5e\n', n, seqtime, partime);
  end

  function res = parfaster_tictoc(n, fpopts, ntests, nthreads, fpclass)
    X = rand(n, 1, fpclass);
    seqtimings = zeros(1, ntests);
    partimings = zeros(1, ntests);
    for i = 1:ntests
      tic;
      Y = cpfloat(X, fpopts, 1);
      seqtimings(i) = toc();
      tic;
      Y = cpfloat(X, fpopts, -nthreads);
      partimings(i) = toc();
    end
    seqtime = median(seqtimings);
    partime = median(partimings);
    res = partime < seqtime;
    fprintf('[%7d]   %.5e   %.5e\n', n, seqtime, partime);
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
