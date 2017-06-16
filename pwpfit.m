function [fitobject, x0] = pwpfit (x1, x2, y, n, x0)
%PWPFIT Fits piece-wise polynomial functions to data under constraints.
%
% Finds a piece-wise defined, polynomial function
%
%   f(x) = f1(x) if x <= x0, f2(x) else,
%
% where f1, f2 are polynomials in x of degree n; 
% minimizing
%
%   sum[j=1:k1] |f1(x1(j)) - y(j)|^2 + sum[j=1:k2] |f2(x2(j)) - y(k1+j)|^2,
%
% where k1, k2 are the length of x1, x2, respectively, and k1+k2 = n is the
% length of y;
% subject to
%
%   f1(x0) == f2(x0).
%
%% Usage and description
%
%   fitobject, x0 = pwpfit(x1, x2, y, n)
%       (...)     = pwpfit(..., x0)
%
% If there is no |x0| given, it is calculated based on the fit of |f1| and
% |f2|.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-22
% * Changed:    2017-02-23
%
%%

% column of monomials to degree n
% p(x) = [1,...,x^n]^T
p = monomials(n);

% length of piece-wise data
% k1 = #x1 = #y1
k1 = length(x1);
% k2 = #x2 = #y2
k2 = length(x2);

%% Reduction to least-square optimization
%
% As f1, f2 are polynomials of degree n, i.e.
%
%   fi = qin x^n + ... + qi1 x + qi0,           i=1,2,
%
% the objective can be written as least-square problem in q = [q1 q2]^T:
%
%   find q minimizing || C*q - y ||^2,
%
% where ||.|| is the L2-norm and
%
%       | 1 x1,1  ... x1,1^n  |                     |
%   	| :   :    \     :    |          0          |
%       | 1 x1,k1 ... x1,k1^n |                     |
%   C = | ------------------------------------------|,
%       |                     | 1 x2,1  ... x2,1^n  |
%   	|          0          | :   :    \     :    |
%       |                     | 1 x2,k2 ... x2,k2^n |
%
% subject to the equality constraint which is equivalent to the matrix 
% equality
%
%       [1 x0 ... x0^n]*q1 = [1 x0 ... x0^n]*q2
%   <=>
%       [1 -1 x0 -x0 ... x0^n -x0^n]*q = 0.
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
C = zeros(k1+k2, 2*(n+1));
d = y;
for j = 1:k1
    C(j,1+(0:n)) = double(p(x1(j))');
end
for j = 1:k2
    C(k1+j, n+1+(0:n)+1) = double(p(x2(j))');
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
    x0 = fsolve(@(x) double(f1(x)-f2(x)), x1(end));
end

fitobject = pwfitobject(sprintf('poly%g', n), {f1, f2}, x0, [q1 q2], n);

% piece-wise function f
f(x) = piecewise(x<=x0, f1(x), f2(x));

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