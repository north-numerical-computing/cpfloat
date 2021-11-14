% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

rng(1)
sizes = [50:50:1000];
nsizes = length(sizes);

options.precision = 'h';
options.subnormal = 1;
options.round = 1;
chop(1, options);
cpfloat(1, options);

times_chop = zeros(1, nsizes);
times_cpfloat = zeros(1, nsizes);
times_double = zeros(1, nsizes);

for k = 1:nsizes
  n = sizes(k);
  A = cpfloat(randn(n, n), options);

  lu_fun = @()lu_chop(A, @(x)chop(x, options));
  times_chop(k) = timeit(lu_fun);

  lu_fun = @()lu_chop(A, @(x)cpfloat(x, options));
  times_cpfloat(k) = timeit(lu_fun);

  lu_fun = @()lu_chop(A, @(x)(x));
  times_double(k) = timeit(lu_fun);
end

semilogy(sizes, [times_chop; times_cpfloat; times_double]);

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
