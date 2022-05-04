% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

warning('off', 'all')

% Parameters
fasttiming = false;
roundingmodes = [1:4];
intlab_roundingmodes = [0 1 -1 2];
sizes = 1:9;
sizes = [sizes, 10*sizes, 100*sizes, 1000*sizes, 10000];

subnormal = 1;
explim = 0;
flip = 0;

nsizes = length(sizes);
formats = {'h', 'b', 'c'};
nroundingmodes = length(roundingmodes);
nformats = length(formats);
speedups = zeros(nroundingmodes,nsizes,nformats);

for k = 1:nformats
  format = formats{k};
  fprintf("Target format: %s\n", format);

  switch format
    case 'h'
      precision = 11;
      emax = 15;
    case 'b'
      precision = 8;
      emax = 127;
    case 's'
      precision = 23;
      emax = 127;
    case 'd'
      precision = 64;
      emax = 1023;
    case 'c'
      % TensorFloat-32
      precision = 11;
      emax = 127;
  end
  emin = 1 - emax;
  emins = emin + 1 - precision;
  xmins = 2^emins;
  xmin = 2^emin;

  for j = 1:nsizes
    n = sizes(j);
    m = n;
    fprintf("   n = %5d\n", n);

    % Generate test matrix
    A = rand(m, n);
    if generatesubnormals
      A = A * (xmin-xmins) + xmins;
    else
      A = A + xmin;
    end

    Ac = zeros(m, n);
    Acf = zeros(m, n);

    for i = 1:nroundingmodes
      round = roundingmodes(i);
      setround(intlab_roundingmodes(i));

      clear options
      options.format = format;
      options.subnormal = subnormal;
      options.round = round;
      options.flip = flip;
      options.explim = explim;
      if format == 'c'
        options.params = [precision, emax];
      end

      if fasttiming
        t1 = tic;
        Y = cpfloat(A, options);
        time_cpfloat = toc(t1);

        t2 = tic;
        Z = flround(A, precision, emax);
        time_intlab = toc(t2);
      else
        cpfloat_fun = @()cpfloat(A, options);
        time_cpfloat = timeit(cpfloat_fun);

        intlab_fun = @()flround(A, precision, emax);
        time_intlab = timeit(intlab_fun);
      end

      speedups(i,j,k) = time_intlab/time_cpfloat;

    end
  end
end

if generatesubnormals
  fprefix = sprintf('%s/speedup-intlab-subnormals', datdir);
else
  fprefix = sprintf('%s/speedup-intlab-normals', datdir);
end
for k = 1:nformats
  fileid = fopen(sprintf('%s-%s-%c.dat', fprefix, 'double', formats{k}), 'w');
  for j = 1:nsizes
    fprintf(fileid, '%5d ', sizes(j));
    for i = 1:nroundingmodes
      fprintf(fileid, '%20.15e ', speedups(i,j,k));
    end
    fprintf(fileid, '\n');
  end
end

plotcols = 2;
plotrows = ceil(nformats/plotcols);
for i = 1:nformats
  figure(1)
  subplot(plotrows, plotcols, i);
  loglog(sizes,speedups(:,:,i));
  axis([min(sizes),max(sizes),0.1,50]);

  figure(2)
  subplot(plotrows, plotcols, i);
  semilogx(sizes,speedups(:,:,i));
  axis([min(sizes),max(sizes),1,50]);
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
