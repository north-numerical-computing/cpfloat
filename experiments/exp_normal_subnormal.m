% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

warning('off', 'all')

% Parameters
fasttiming = false;
roundingmodes = [-1:7];
sizes = 1:9;
sizes = [10*sizes , 100*sizes, 1000*sizes, 10000];

subnormal = 1;
explim = 0;
flip = 0;

nsizes = length(sizes);
if strcmp(storageformat, 'single')
  formats = {'h', 'b', 'c'};
else
  formats = {'h', 'b', 'c'};
end
nroundingmodes = length(roundingmodes);
nformats = length(formats);
speedups = zeros(nroundingmodes,nsizes,nformats);

for k = 1:nformats
  format = formats{k}

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
      if strcmp(storageformat,'single')
        precision = 11;
        emax = 127;
      else
        precision = 22;
        emax = 1023;
      end
      params = [precision, emax];
      options.params = params;
  end
  emin = 1-emax;
  emins = emin + 1 - precision;
  xmins = 2^emins;
  xmin = 2^emin;

  for j = 1:nsizes
    n = sizes(j);
    m = n;

    % Generate test matrix
    Anormal = rand(m, n, storageformat) + xmin;
    Asubnormal = rand(m, n, storageformat) * (xmin-xmins) + xmins;

    for i = 1:nroundingmodes
      round = roundingmodes(i);

      options.format = format;
      options.subnormal = subnormal;
      options.round = round;
      options.flip = flip;
      options.explim = explim;

      if fasttiming
        t1 = tic;
        Y = cpfloat(Anormal, options);
        time_normal= toc(t1);

        t2 = tic;
        Z = cpfloat(Asubnormal, options);
        time_subnormal= toc(t2);
      else
        cpfloat_fun = @()cpfloat(Anormal, options);
        time_normal= timeit(cpfloat_fun);

        chop_fun = @()cpfloat(Asubnormal, options);
        time_subnormal= timeit(chop_fun);
      end

      speedups(i,j,k) = time_subnormal/time_normal;

    end
  end
end

fprefix = sprintf('%s/speedup-normals-subnormals', datdir);

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

plotcols = 3;
plotrows = ceil(nformats/plotcols);
for i = 1:nformats
  figure(2)
  subplot(plotrows, plotcols, i);
  semilogx(sizes,speedups(:,:,i));
  axis([min(sizes),max(sizes),0,2]);
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
