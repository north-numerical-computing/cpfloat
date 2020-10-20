% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

if (~exist(sprintf('%s/overhead-clang-single.dat', datdir), 'file') ...
    || ~exist(sprintf('%s/overhead-clang-double.dat', datdir), 'file'))
  system('make -C .. run_exp_overhead');
end

% Input parameters.
sizes = 1:9;
sizes = [sizes*100 sizes*1000 10000];
nsizes = length(sizes);
roundingmodes = [-1:7];
nroundings = length(roundingmodes);

ntests = 10;

options.format = "h";
options.subnormal = 1;
options.flip = 0;
options.explim = 1;
fmin = 2^-14; % smallest positive normal in binary16

storageformats = {'single', 'double'};
fignumber = 0;
for iter1 = 1:2

  storageformat = storageformats{iter1};
  matlabdata = sprintf('%s/overhead-matlab-%s.dat', datdir, storageformat);

  if (~exist(matlabdata, 'file'))

    fid = fopen(matlabdata, 'w');

    medtimings = zeros(nsizes,nroundings);
    timing = zeros(1,ntests);
    for i = 1:nsizes
      n = sizes(i);
      X = rand(n,n,storageformat) + cast(fmin,storageformat);
      assert(all(class(X) == storageformat))
      fprintf(fid, "%5d", n);
      fprintf("%5d", n);
      for round = 1:nroundings
        options.round = roundingmodes(round);
        for k = 1:ntests
          start = tic;
          Y = cpfloat(X, options);
          timing(k) = toc(start);
        end
        timing = sort(timing);
        medtimings(i, round) = timing(floor(ntests/2));
        fprintf(fid, " %10.5e", medtimings(i, round));
        fprintf(" %10.5e", medtimings(i, round));
      end
      fprintf(fid, "\n");
      fprintf("\n");
    end

    fclose(fid);
  end

  cdata = sprintf('%s/overhead-clang-%s.dat', datdir, storageformat);

  ctimings = textread(cdata);
  matlabtimings = textread(matlabdata);
  assert(all(ctimings(:,1) == matlabtimings(:,1)));

  speedups = matlabtimings(:,2:end) ./ ctimings(:,2:end);

  resfile = sprintf('%s/overhead-%s.dat', datdir, storageformat);
  fires = fopen(resfile, 'w');
  for i = 1:nsizes
    n = sizes(i);
    fprintf(fires, "%5d", n);
    for j = 1:nroundings
      fprintf(fires, " %10.5e", speedups(i, j));
    end
    fprintf(fires, "\n");
  end
  fclose(fires);

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
