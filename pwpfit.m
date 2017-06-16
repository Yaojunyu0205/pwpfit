function [fitobject, x0] = pwpfit (xa, xb, z, n, x0)
%PWPFIT Fits piece-wise polynomial functions to data under constraints.
%
% Finds a piece-wise defined, polynomial function
%
%   f(x) = fa(x1,...,xm) if x1 <= x0, fb(x1,...,xm) else,
%
% where fa, fb are polynomials in x1,...,xm of degree n,
%
%   fa(x) = an0 x1^n + ... + a0n xm^n + ... + a10 x1 + ... + a01 xm + a0;
%   fb(x) = bn0 x1^n + ... + b0n xm^n + ... + b10 x1 + ... + b01 xm + b0;
%
% minimizing
%
%   sum[j=1:ka] |fa(xa1(j),...,xam(j)) - y(j)|^2 
%                        + sum[j=1:kb] |fb(xb1(j),...,xbm(j)) - y(ka+j)|^2,
%
% where ka, kb are the length of xa, xb, respectively, and ka+kb = k is the
% length of y;
% subject to
%
%   fa(x0,...) == fb(x0,...)
%
% for all x2,...,xm in R^(m-1).
%
%% Usage and description
%
%   [fitobject, x0] = pwpfit(xa, xb, y, n)
%   [...] = pwpfit(..., x0)
%
% Returns fit of xa, xb against y, where xa, xb, xy are column vectors with
% size([xa; xb]) = size(y).
% If there is no |x0| given, it is calculated based on the fit of fa and
% fb.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-22
% * Changed:    2017-06-16
%
%%

assure(size(xa, 2) == size(xb, 2), 'xa and xb must have same number of columns.');
assure(all(size([xa; xb]) == size(z)), '[xa; xb] and y must equal in size.');

% number of columns in data
m = size(xa, 2);

% column of monomials to degree n
% p(x) = [1,...,x^n]^T
[p, X, r] = monomials(n, m);

% length of piece-wise data
% k1 = #x1 = #y1
ka = length(xa);
% k2 = #x2 = #y2
kb = length(xb);


%% Reduction to least-square optimization
%
% As fa, fb are polynomials of degree n, i.e.
%
%   fi = qin0 x1^n + ... + qi0n xm^n + ... + qi10 x1 + ... + qi01 xm + qi0,
%
% the objective can be written as least-square problem in q = [q1 q2]^T:
%
%   find q minimizing || C*q - y ||^2,
%
% where ||.|| is the L2-norm and
%
%       | Ca |    |
%   C = |---------|
%       |    | Cb |
%
% with
%
%        | 1 xi1,1  ... xim,1  ... xi1,1^n  ... xim,1^n   |
%   Ci = | :    :    \     :    \     :      \     :      |
%        | 1 xi1,ki ... xim,ki ... xi1,ki^n ... xim,ki^n  |
%
% subject to the curve equality constraint fa(x0) == fb(x0), 
% which for m = 1 is equivalent to the matrix equality
%
%       [1 x0 ... x0^n]*q1 = [1 x0 ... x0^n]*q2
%   <=>
%       [1 -1 x0 -x0 ... x0^n -x0^n]*q = 0.
%
% For m > 1, we have the surface equality constraint
%
%   fa(x0, X') = fb(x0, X')
%
% for all X' in R^(m-1) equivalent to the matrix equality
%
%       []
%

if nargin < 5
    % no equality constraint
    Aeq = [];
    beq = [];
    x0 = NaN;
else
    % equality constraint
    % Aeq1*q1 + Aeq2*q2 = beq
    Aeq1 = double(p(x0)');
    Aeq = [Aeq1 -Aeq1];
    beq = 0;
end

% least squares objective
% find q minimizing the L2-norm
% ||C*q-d||^2
C = zeros(ka+kb, 2*(n+1));
d = z;
for j = 1:ka
    C(j,1+(0:n)) = double(p(xa(j))');
end
for j = 1:kb
    C(ka+j, n+1+(0:n)+1) = double(p(xb(j))');
end

% inequality condition
% A*q <= b
A = ones(1,2*(n+1));
b = 1e4;

% solve LSQ for q
q = lsqlin(C, d, A, b, Aeq, beq);

% piece-wise coefficients
q1 = q(1+(0:n));
q2 = q(n+1+(0:n)+1);

% piece-wise functions
syms x
f1(x) = q1'*p(x);
f2(x) = q2'*p(x);

% if no x0 was given, find x0 s.t. f1(x0) == f2(x0)
if isnan(x0)
    x0 = fsolve(@(x) double(f1(x)-f2(x)), xa(end));
end

fitobject = pwfitobject(sprintf('poly%g', n), {f1, f2}, x0, [q1 q2], n);

% piece-wise function f
% f(x) = piecewise(x<=x0, f1(x), f2(x));

end

% function p = monomials(n)
% %MONOMIALS Creates a column vector of monomials to degree n.
% %   Vector p is symbolic function of x, i.e. p(x) = [1,...,x^n]^T.
%     syms x
%     P = sym('P', [n 1]);
%     for i = 0:n
%         P(1+i) = x^i;
%     end
%     p = symfun(P, x);
% end