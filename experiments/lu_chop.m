% SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
% SPDX-License-Identifier: LGPL-2.1-or-later

function [L,U,p] = lu_chop(A, chop_fun)
%LUTX_CHOP LU factorization with partial pivoting.
%    [L,U,P] = LU_CHOP(A,CHOP_FUN) computes the LU factorization of an N-by-N
%    matrix A using Gaussian Elimination with partial pivoting, and CHOP_FUN to
%    round the result of each operation performed during the factorization.
%    The matrix L and U are unit lower triangular and upper triangular,
%    respectively, and P is a permutation of the vector 1:N, which satisfy
%    L * U = A(P,:).

  if nargin < 1 || nargin > 2
    error('lu_chop:invalidNargin', 'This functions accepts one or two input arguments.');
  end
  if nargin < 2
    chop_fun = @(x) cpfloat(x);
  else
    if (~isa(chop_fun, 'function_handle'))
      error('lu_chop:invalidRoundingFunction',...
            'The second argument must be a function handle');
    end
  end

  [m,n] = size(A);
  if (~isnumeric(A) || m ~= n)
    error('lu_chop:invalidMatrixSize',...
          'The first argument must be a square matrix.');
  end
  p = [1:n];

  for k = 1:n-1

    % Find pivot.
    [x,j] = max(abs(A(k:n, k)));
    j = j+k-1;
    if (x < eps)
      error('lu_chop:singular', 'Breakdown during Gaussian elimination.');
    else
      % Swaps rows.
      if (j ~= k)
        A([j k],:) = A([k j], :);
        p([j k]) = p([k j]);
      end

      % Compute column of L.
      A(k+1:n,k) = chop_fun(A(k+1:n,k) / A(k,k));

      % Update trailing block of A.
      A(k+1:n,k+1:n) = chop_fun(A(k+1:n,k+1:n) - chop_fun(A(k+1:n,k) * A(k,k+1:n)));
    end
  end

  % Extract output
  L = tril(A,-1);
  L(1:n+1:n*n) = 1;
  U = triu(A);
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
