% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

function cpfloat_test
%TEST_CPFLOAT Test the cpfloat function.
%   The tests are for single precision and fp16.

  clear cpfloat fp options options2 assert_eq

  usingoctave = exist('OCTAVE_VERSION', 'builtin');

  if usingoctave
    rand('seed', 1);
  else
    rng(1);
  end

  n = 0;
  uh = 2^(-11);  % Unit roundoff for fp16.
  pi_h = 6432*uh; % fp16(pi)

  % Check handling of defaults and persistent variable.
  fp.format = 'bfloat16'; [c,options] = cpfloat(pi,fp);
  assert_eq(fp.format,options.format)
  assert_eq(options.subnormal,0)

  fp.format = []; [c,options] = cpfloat(pi,fp);
  assert_eq(options.format,'h')  % Check default;

  fp.subnormal = 0; [c,options] = cpfloat(pi,fp);
  assert_eq(options.subnormal,0)

  fp.subnormal = []; [c,options] = cpfloat(pi,fp);
  assert_eq(options.subnormal,1)  % Check default;

  fp.round = []; [c,options] = cpfloat(pi,fp);
  assert_eq(options.round,1)  % Check default.

  fp.flip = []; [c,options] = cpfloat(pi,fp);
  assert_eq(options.flip,0)  % Check no default.

  fp.explim = []; [c,options] = cpfloat(pi,fp);
  assert_eq(options.explim,1)  % Check default.

  fp.explim = 0; [c,options] = cpfloat(pi,fp);
  assert_eq(options.explim,0)  % Check no default.

  clear cpfloat fp options
  fp.flip = 1; [~,options] = cpfloat([],fp);
  assert_eq(options.format,'h')
  assert_eq(options.round,1)
  assert_eq(options.subnormal,1)

  clear cpfloat fp options
  % check all default options
  fp.format = []; fp.subnormal = [];
  fp.round = []; fp.flip = [];
  fp.p = [];
  [c,options] = cpfloat(pi,fp);
  assert_eq(options.format,'h')
  assert_eq(options.subnormal,1)
  assert_eq(options.round,1)
  assert_eq(options.flip,0)
  assert_eq(options.p,0.5)
  % % Takes different path from previous test since fpopts exists.
  % fp.subnormal = 0;
  % fp.format = []; [c,options] = cpfloat(pi,fp);
  % assert_eq(options.format,'h')

  % Check flip output.
  clear cpfloat fp
  fp.flip = 1; fp.format = 'd';
  c = ones(8,1);
  d = cpfloat(c,fp); assert_eq(norm(d-c,1)>0,true);
  d = cpfloat(c',fp); assert_eq(norm(d-c',1)>0,true);
  fp.p = 0; % No bits flipped.
  d = cpfloat(c,fp); assert_eq(d,d);
  fp.p = 1; % All bits flipped.
  d = cpfloat(c,fp); assert_eq(all(d ~= c),true);

  clear cpfloat
  [~,fp] = cpfloat;
  assert_eq(fp.subnormal,1)
  assert_eq(fp.format,'h')
  [c,options] = cpfloat(pi);
  assert_eq(options.format,'h')
  assert_eq(options.subnormal,1)
  assert_eq(options.round,1)
  assert_eq(options.flip,0)
  assert_eq(options.p,0.5)

  clear fp
  fp.format = 'd';
  [c,options] = cpfloat(pi,fp);
  assert_eq(options.format,'d')
  assert_eq(options.subnormal,1)
  assert_eq(options.params, [53 1023])
  [~,fp] = cpfloat;
  assert_eq(fp.format,'d')
  assert_eq(fp.subnormal,1)
  assert_eq(fp.params, [53 1023])

  clear fp
  fp.format = 'bfloat16'; [c,options] = cpfloat(pi,fp);
  assert_eq(options.format,'bfloat16')
  assert_eq(options.subnormal,0)
  assert_eq(options.params, [8 127])
  [~,fp] = cpfloat;
  assert_eq(fp.format,'bfloat16')
  assert_eq(fp.subnormal,0)
  assert_eq(fp.params, [8 127])

  clear cpfloat
  [~,fp] = cpfloat;
  fp.format = 'b';
  fp = rmfield(fp, 'params');
  [c,options] = cpfloat(pi,fp);
  assert_eq(options.subnormal,1) % No subnormals only if that field was empty.

  % Check these usages do not give an error.
  c = cpfloat([]);
  cpfloat([]);
  cpfloat([],fp);
  cpfloat(1,[]);
  cpfloat(1,fp);
  c = cpfloat(1,fp);

  % Test matrix.
  options.format = 'b';
  options = rmfield(options, 'params');
  A = magic(4);
  C = cpfloat(A,options);
  assert_eq(A,C);
  B = A + randn(size(A))*1e-12;
  C = cpfloat(B,options);
  assert_eq(A,C);
  A2 = hilb(6); C = cpfloat(A2);

  options.format = 'c';
  options.params = [8 127];  % bfloat16
  C1 = cpfloat(A,options);
  assert_eq(A,C1);
  C2 = cpfloat(B,options);
  assert_eq(A,C2);
  assert_eq(C,cpfloat(A2));

  clear options
  options.format = 'c';
  options.params = [11 15];  % h
  options2.format = 'h';
  A = hilb(6);
  [X1,opt] = cpfloat(A,options);
  [X2,opt2] = cpfloat(A,options2);
  assert_eq(X1,X2)
  % assert_eq(cpfloat(A,options),cpfloat(A,options2));

  % Row vector
  clear options
  options.format = 'h';
  A = -10:10;
  C = cpfloat(A,options);
  assert_eq(A,C);
  B = A + randn(size(A))*1e-12;
  C = cpfloat(B,options);
  assert_eq(A,C);

  % Column vector
  options.format = 's';
  A = (-10:10)';
  C = cpfloat(A,options);
  assert_eq(A,C);
  B = A + A.*rand(size(A))*1e-14;  % Keep 0 as 0.
  C = cpfloat(B,options);
  assert_eq(A,C);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Main loop: test single and half formats.
  for i = 1:4
    clear cpfloat fp options

    if i == 1
      % Single precision tests.
      [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('single');
      options.format = 's';
    elseif i == 2
      % Half precision tests.
      [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('half');
      options.format = 'h'
    elseif i == 3
      % Quarter precision tests.
      [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('q43');
      options.format = 'E4M3';
    elseif i == 4
      % Quarter precision tests.
      [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('q52');
      options.format = 'E5M2';
    end
    options.subnormal = 0;

    x = pi;
    if i == 1
      y = double(single(x));
    elseif i == 2
      y = pi_h; % double(fp16(x));
    elseif i == 3
      y = 3.25;
    elseif i == 4
      y = 3.0;
    end
    c = cpfloat(x,options);
    assert_eq(c,y);
    x = -pi;
    c = cpfloat(x,options);
    assert_eq(c,-y);

    % Next number power of 2.
    y = 2^10;
    if i == 1
      dy = double(eps(single(y)));
    elseif i == 2
      dy = 2*y*uh; % double(eps(fp16(y)));
    elseif i == 3
      y = 2^4;
      dy = 2*y*u;
    elseif i == 4
      y = 2^4;
      dy = 2*y*u;
    end
    x = y + dy;
    c = cpfloat(x,options);
    assert_eq(c,x)

    % Number just before a power of 2.
    x = y - dy;
    c = cpfloat(x,options);
    assert_eq(c,x)

    % Next number power of 2.
    y = 2^(-4);
    if i == 1
      dy = double(eps(single(y)));
    elseif i == 2
      dy = 2*y*uh; % double(eps(fp16(y)));
    elseif i == 3
      dy = 2*y*u;
    elseif i == 4
      dy = 2*y*u;
    end
    x = y + dy;
    c = cpfloat(x,options);
    assert_eq(c,x)

    % Check other rounding options
    for rmode = 1:6
      options.round = rmode;
      x = y + (dy*10^(-3));
      c = cpfloat(x,options);
      if options.round == 2
        assert_eq(c,y+dy) % Rounding up.
      elseif options.round >= 5
        % Check rounded either up or down.
        if c ~= y+dy
          assert_eq(c,y);
        end
      else
        assert_eq(c,y);
      end
    end

    % Overflow tests.
    for j = 1:6
      options.round = j;
      x = xmax;
      c = cpfloat(x,options);
      assert_eq(c,x)
    end

    % Infinities tests.
    for j = 1:6
      options.round = j;
      x = inf;
      c = cpfloat(x,options);
      assert_eq(c,x)
      c = cpfloat(-x,options);
      assert_eq(c,-x)
    end

    % IEEE 754-2019, page 27: rule for rounding to infinity.
    % Round to nearest
    options.round = 1; % reset the rounding mode to default
    x = 2^emax * (2-(1/2)*2^(1-p));  % Round to inf.
    c = cpfloat(x,options);
    assert_eq(c,inf)
    c = cpfloat(-x,options);
    assert_eq(c,-inf)

    x = 2^emax * (2-(3/4)*2^(1-p));  % Round to realmax.
    c = cpfloat(x,options);
    assert_eq(c,xmax)
    c = cpfloat(-x,options);
    assert_eq(c,-xmax)

    % Round toward plus infinity
    options.round = 2;
    x = 2^emax * (2-(1/2)*2^(1-p));
    c = cpfloat(x,options);
    assert_eq(c,inf)
    c = cpfloat(-x,options);
    assert_eq(c,-xmax)

    % Round toward minus infinity
    options.round = 3;
    c = cpfloat(x,options);
    assert_eq(c,xmax)
    c = cpfloat(-x,options);
    assert_eq(c,-inf)

    % Round toward zero
    options.round = 4;
    c = cpfloat(x,options);
    assert_eq(c,xmax)
    c = cpfloat(-x,options);
    assert_eq(c,-xmax)

    % Round to nearest.
    options.round = 1; % reset the rounding mode to default
    if i == 2
      x = 1 + 2^(-11);
      c = cpfloat(x,options);
      assert_eq(c,1)
    end

    % Underflow tests.
    if i == 1
      delta = double(eps(single(1)));
    else
      delta = 2*uh; % double(eps(fp16(1)));
    end

    options.subnormal = 1;
    c = cpfloat(xmin,options); assert_eq(c,xmin)
    x = [xmins xmin/2 xmin 0 xmax 2*xmax 1-delta/5 1+delta/4];
    c = cpfloat(x,options);
    c_expected = [x(1:5) inf 1 1];
    assert_eq(c,c_expected)

    options.subnormal = 0;
    c = cpfloat(xmin,options); assert_eq(c,xmin)
    x = [xmins xmin/2 xmin 0 xmax 2*xmax 1-delta/5 1+delta/4];
    c = cpfloat(x,options);
    c_expected = [0 0 x(3:5) inf 1 1];
    assert_eq(c,c_expected)

    % Smallest normal number and spacing between the subnormal numbers.
    y = xmin; delta = xmin*2^(1-p);
    x = y - delta; % The largest subnormal number.
    options.subnormal = 1;
    c = cpfloat(x,options);
    assert_eq(c,x)
    % Round up if subnormals are not supported.
    options.subnormal = 0;
    c = cpfloat(x,options);
    assert_eq(c,xmin)
    % Flush subnormals to zero if subnormals are not supported.
    options.subnormal = 0;
    c = cpfloat(xmins,options);
    assert_eq(c,0)

    options.subnormal = 1;
    x = xmins*8;  % A subnormal number.
    c = cpfloat(x,options);
    assert_eq(c,x)

    % Numbers smaller than smallest representable number.
    options.subnormal = 0;
    x = xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,-0)
    x = xmin / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmin / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)

    options.subnormal = 1;
    x = xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)

    % Do not limit exponent.
    options.explim = 0;
    x = xmin/2;  c = cpfloat(x,options); assert_eq(c,x)
    x = -xmin/2;  c = cpfloat(x,options); assert_eq(c,x)
    x = xmax*2;  c = cpfloat(x,options); assert_eq(c,x)
    x = -xmax*2;  c = cpfloat(x,options); assert_eq(c,x)
    x = xmins/2; c = cpfloat(x,options); assert_eq(c,x)
    x = -xmins/2; c = cpfloat(x,options); assert_eq(c,x)
    A = [pi -pi; pi -pi];
    C = cpfloat(A,options);
    options.explim = 1;
    assert_eq(C,cpfloat(A,options));

    % Round toward plus infinity
    options.round = 2;
    options.subnormal = 0;
    x = xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,xmin)
    x = -xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)

    options.subnormal = 1;
    x = xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,xmins)
    x = -xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,xmins)
    x = -xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)

    % Round toward minus infinity
    options.round = 3;
    options.subnormal = 0;
    x = xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,-xmin)

    options.subnormal = 1;
    x = xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,-xmins)
    x = xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,-xmins)

    % Round toward zero.
    options.round = 4;
    options.subnormal = 0;
    x = xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmin / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)

    options.subnormal = 1;
    x = xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmins / 2;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)
    x = -xmins / 4;
    c = cpfloat(x,options);
    assert_eq(c,0)

  end % for i
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clear options

  % Test rounding with CHOPFAST versus native rounding.
  options.format = 's';
  m = 100; y = zeros(3,n); z = y;
  for i = 1:m
    x = randn;
    options.round = 2; y(i,1) = cpfloat(x,options);
    options.round = 3; y(i,2) = cpfloat(x,options);
    options.round = 4; y(i,3) = cpfloat(x,options);
    if usingoctave
      fesetround(inf); z(i,1) = single(x);
      fesetround(-inf); z(i,2) = single(x);
      fesetround(0); z(i,3) = single(x);
    else
      % Use undocumented function to set rounding mode in MATLAB.
      feature('setround',inf), z(i,1) = single(x);
      feature('setround',-inf), z(i,2) = single(x);
      feature('setround',0), z(i,3) = single(x);
    end
  end
  assert_eq(y,z)
  % Switch back to round to nearest.
  if usingoctave
    fesetround(0.5);
  else
    feature('setround',0.5)
  end

  % Double precision tests.
  [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('d');
  options.format = 'd';
  x = [1e-309 1e-320 1 1e306];  % First two entries are subnormal.
  c = cpfloat(x,options);
  assert_eq(c,x)
  options.subnormal = 0;
  c = cpfloat(x,options);
  assert_eq(c,[0 0 x(3:4)])

  options.format = 'd'; options.subnormal = 0; cpfloat([],options);
  a = cpfloat(pi); assert_eq(a,pi)
  options.format = 'd'; options.subnormal = 1; cpfloat([],options);
  a = cpfloat(pi); assert_eq(a,pi)

  x = pi^2;
  clear options
  options.format = 'd';
  y = cpfloat(x,options);  % Should not change x.
  assert_eq(x,y);
  options.round = 2;
  y = cpfloat(x,options);  % Should not change x.
  assert_eq(x,y);
  options.round = 3;
  y = cpfloat(x,options);  % Should not change x.
  assert_eq(x,y);
  options.round = 4;
  y = cpfloat(x,options);  % Should not change x.
  assert_eq(x,y);

  % Test on single inputs.
  clear options
  ps = single(pi);
  pd = double(ps);
  options.format = 'b';
  ys = cpfloat(ps,options);
  assert_eq(isa(ys,'single'),true)
  yd = cpfloat(pd);
  assert_eq(double(ys),yd)

  options.format = 'h'; options.round = 2;
  as = single(rand(n,1)); ad = double(as);
  delta = single(rand(n,1));
  cd = cpfloat(ad + 1e-5*double(delta),options);
  cs = cpfloat(as + 1e-5*delta,options);
  assert_eq(cd,double(cs));

  options.format = 'c';
  options.params = [11 5];
  temp1 = cpfloat(single(pi),options);
  options.format = 'h';
  options = rmfield(options, 'params');
  temp2 = cpfloat(single(pi),options);
  assert_eq(temp1,temp2)

  % Test base 2 logarithm
  options.format = 'h';
  options.round = 4;
  x = single(2^-3 * (sum(2.^(-[0:23]))));
  assert_eq(cpfloat(x,options), single(2^-3 * (sum(2.^(-[0:10])))))

  x = 2^-3 * (sum(2.^(-[0:52])));
  assert_eq(cpfloat(x,options), 2^-3 * (sum(2.^(-[0:10]))))

  options.format = 's';
  x = single(2^-3 * (sum(2.^(-[0:23]))));
  assert_eq(cpfloat(x,options), x)

  x = 2^-3 * (sum(2.^(-[0:52])));
  assert_eq(cpfloat(x,options), 2^-3 * (sum(2.^(-[0:23]))))

  options.format = 'd';
  x = 2^-3 * (sum(2.^(-[0:52])));
  assert_eq(cpfloat(x,options), x)

  options.round = 1;
  temp = 0;
  try
    options.format = 'c';
    options.params = [12 5];
    temp = cpfloat(single(pi),options); % Error - double rounding!
  catch
  end
  assert_eq(temp,0)
  try
    options.format = 'c';
    options.params = [26 9];
    temp = cpfloat(pi,options); % Error - double rounding!
  catch
  end
  assert_eq(temp,0)
  try
    temp = cpfloat(complex(1,1)); % Error - complex data!
  catch
  end
  assert_eq(temp,0)

  fprintf('All tests successful!\n')

  clear cpfloat fp options options2 assert_eq

  %%%%%%%%%%%%%%%%%%%%%%%
  function assert_eq(a,b)
  % if isempty(n), n = 0; end  % First call.
    n = n+1;
    if ~isequal(a,b)
      error('Failure')
    end
    fprintf('Test %g succeeded.\n', n )
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
