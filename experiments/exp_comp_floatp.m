% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

warning('off', 'all')

% Parameters
fasttiming = false;
roundingmodes = [1:6];
% storageformat= 'single';
storageformat= 'double';
% generatesubnormals = true;
sizes = [1:1:10, 20:10:100];

subnormal = 1;
explim = 0;
flip = 0;

nsizes = length(sizes);
if strcmp(storageformat, 'single')
  formats = {'h', 'b', 'c'};
else
  formats = {'h', 'b', 'c'};
end
% Floatp package requires numbers of bits specified
exponents = [5, 8, 8];
significands = [10, 7, 10];
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
      % Largest format allowed
      % if strcmp(storageformat,'single')
      %   precision = 11;
      %   emax = 127;
      % else
      %   precision = 22;
      %   emax = 1023;
      % end
      % TensorFloar-32
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
    A = rand(m, n, storageformat);
    if generatesubnormals
      A = A * (xmin-xmins) + xmins;
    else
      A = A + xmin;
    end

    Ac = zeros(m, n, storageformat);
    Acf = zeros(m, n, storageformat);

    for i = 1:nroundingmodes
      round = roundingmodes(i);

      clear options
      options.format = format;
      options.subnormal = subnormal;
      options.round = round;
      options.flip = flip;
      options.explim = explim;
      if format == 'c'
        options.params = [precision, emax];
      end

      f_d_init_bits_expo(exponents(k));
      f_d_init_round(round);

      if fasttiming
        t1 = tic;
        Y = cpfloat(A, options);
        time_cpfloat = toc(t1);

        t2 = tic;
        Z = f_d_dec2floatp(A, significands(k));
        time_floatp = toc(t2);
      else
        cpfloat_fun = @()cpfloat(A, options);
        time_cpfloat = timeit(cpfloat_fun);

        floatp_fun = @()f_d_dec2floatp(A, significands(k));
        time_floatp = timeit(floatp_fun);
      end

      speedups(i,j,k) = time_floatp/time_cpfloat;

    end
  end
end

if generatesubnormals
  fprefix = sprintf('%s/speedup-floatp-subnormals', datdir);
else
  fprefix = sprintf('%s/speedup-floatp-normals', datdir);
end
for k = 1:nformats
  fileid = fopen(sprintf('%s-%s-%c.dat', fprefix, storageformat, formats{k}), 'w');
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
  axis([min(sizes),max(sizes),1,1000]);

  figure(2)
  subplot(plotrows, plotcols, i);
  semilogx(sizes,speedups(:,:,i));
  axis([min(sizes),max(sizes),1,1000]);
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
