function Y = euler(varargin)
% EULER  Euler numbers and Euler polynomials.
%    Y = EULER(N, X) is the N-th Euler polynomial.
%    It is defined as
%    euler(n, x) = subs(diff(2*exp(x*t)/(exp(t) + 1), t, n), t, 0)
%    Y = EULER(N) is the N-th Euler number.
%    It is defined as
%    EULER(n) = 2^n*euler(n,1/2).
%    See also BERNOULLI. 

%   Copyright 2013 The MathWorks, Inc.

narginchk(1, 2);
Y = sym.useSymForNumeric(@euler, varargin{:});
    
