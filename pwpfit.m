function [fitobject, x0] = pwpfit (xa, xb, z, n, x0, varargin)
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

assert(size(xa, 2) == size(xb, 2), 'xa and xb must have same number of columns.');
assert(size([xa; xb],1) == size(z,1), '[xa; xb] and y must have same number of rows.');

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

% zero equality constraint
if ~isempty(varargin) && isnumeric(varargin{1})
    y0 = varargin{1};
    varargin(1) = [];
else
    y0 = NaN;
end



%% Reduction to least-square optimization
%
% As fa, fb are polynomials of degree n, i.e.
%
%   fi = qi0 + qi10 x1 + ... + qi01 xm + ... + qin0 x1^n + ... + qi0n xm^n,
%
% with
%
%   qi = [qi0 qi10 ... qi01 ... qin0 ... qi0n]^T
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
% For m = 2, we have the surface equality constraint
%
%   fa(x0, y) = fb(x0, y)
%
% for all y in R equivalent to the matrix equality
%
%       [1, ..., y^n]*A'*q1 = [y^n, ..., 1]*A'*q2
%   <=>
%       [A' -A']*q = [0, ..., 0]^T
%
% with
%
%        |  1  |  ...  | x0^(n-1) ... 0 | x0^n    0     ... 0 |
%        |-----|       |           \    |      x0^(n-1)       |
%   A' = |-------------|              1 |                \    |.
%        |------------------------------|                   1 |
%        |----------------------------------------------------|
%


%% Curve equality constraint
% Aeq1*q1 - Aeq1*q2 = 0
if ~exist('x0', 'var') || isnan(x0)
    % no equality constraint
    Aeq = [];
    beq = [];
    x0 = NaN;
elseif m == 1
    Aeq1 = double(p(x0)');
    Aeq = [Aeq1 -Aeq1];
    beq = 0;
elseif m == 2
    Aeq1 = zeros(n+1,r);
    j = 0;
    for N=0:n
        [pN, ~, rN] = monomials(N, 1);
        pNx0 = double(pN(x0)');
        Aeq1(1:rN,j+(1:rN)) = diag(pNx0(rN:-1:1));
        j = j + rN;
    end
    Aeq = [Aeq1 -Aeq1];
    beq = zeros(n+1,1);
else
    error('Curve equality constraint for more than 2 variables is not supported yet.');
end

%% Zero constraint
% Aeq*q = 0
if isempty(y0) || all(isnan(y0))
    % no constraint
    Azero = [];
    bzero = [];
elseif m == 1
    Azero1 = double(p(y0)');
    bzero = [0; 0];
    Azero = [Azero1 zeros(1,r); zeros(1,r) Azero1];
elseif m == 2
    Azero1 = zeros(n+1,r);
    j = 0;
    for N=0:n
        [pN, ~, rN] = monomials(N, 1);
        pNy0 = double(pN(y0)');
        Azero1(1:rN,j+(1:rN)) = cdiag(pNy0(rN:-1:1));
        j = j + rN;
    end
    bzero = zeros(2*(n+1),1);
    Azero = [Azero1 zeros(n+1,r); zeros(n+1,r) Azero1];
else
    error('Zero constraint for more than 2 variables is not supported yet.');
end


%% least squares objective
% find q minimizing the L2-norm
% ||C*q-d||^2
C = zeros(ka+kb, 2*r);
d = z;
for j = 1:ka
    Xj = num2cell(xa(j,:));
    C(j,1:r) = double(p(Xj{:})');
end
for j = 1:kb
    Xj = num2cell(xb(j,:));
    C(ka+j, r+(1:r)) = double(p(Xj{:})');
end

% inequality condition
% A*q <= b
A = ones(1,2*r);
b = 1e4;

% solve LSQ for q
q = lsqlin(C, d, A, b, [Aeq; Azero], [beq; bzero]);

% piece-wise coefficients
qa = q(0+(1:r));
qb = q(r+(1:r));

% piece-wise functions
P = formula(p);
Fa = qa'*P;
Fb = qb'*P;
fa = symfun(Fa, X);
fb = symfun(Fb, X);

% if no x0 was given, find x0 s.t. f1(x0) == f2(x0)
if isnan(x0)
    Y0 = num2cell(zeros(1,m-1));
    x0 = fsolve(@(x) double(fa(x,Y0{:})-fb(x,Y0{:})), xa(end));
end

fitobject = pwfitobject(['poly' sprintf('%g', n+zeros(1,m))], {fa, fb}, x0, [qa qb], n, varargin{:});

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