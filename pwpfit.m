function [fitobject, x0, gof, time] = pwpfit (xa, xb, z, n, x0, varargin)
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
%   sum[i=1:ka] |fa(xa1(i),...,xam(i)) - z(i)|^2 
%                        + sum[i=1:kb] |fb(xb1(i),...,xbm(i)) - z(ka+i)|^2,
%
% where ka, kb are the length of xa, xb, respectively, and ka+kb = k is the
% length of z;
% subject to
%
%   fa(x0,...) == fb(x0,...)
%
% for all x2,...,xm in R^(m-1).
%
%% Usage and description
%
%   [fitobject, x0] = pwpfit(xa, xb, z, n)
%   [...] = pwpfit(..., {x0 | NaN}, [y0 | NaN], [pwfoargs...])
%   [..., gof] = pwpfit(...)
%
% Returns fit of xa, xb against z, where xa, xb, xy are column vectors with
% size([xa; xb]) = size(y) and fa(x0,...) == fb(x0,...).
% If there is no x0 given or |x0 == NaN|, x0 is calculated based on the fit
% of fa and fb.
%
% If the optional parameter y0 (and |y0 != NaN|) is given, the returned fit
% is zero in the parameter xj if and only if the j-th component of y0 is
% zero; i.e.
%
%   fa(...,xj=0,...) = fb(...,xj=0,...) = 0
%
% for all x1,...,x[j-1],x[j+1],...,xm in R^(m-1).
%
% The optional argument(s) |pwfoargs| are applied to |fitobject|.
%
% Returns furthermore:
% * |gof| is the goodness-of-fit structure with |gof.rmse| being the Root
% Mean Squared Error:
%
%   rmse = sqrt(sum[i=1:k] |f(x(i)) - z(i)|^2).
%
% * |time| is structure of elapsed execution time where |time.all| is the
% overall execution time; |time.lsq| is time to solve LSQ problem;
% |time.obj| is time to construct LSQ objective matrizes; |time.eq| is time 
% to construct continuity constraint matrix; |time.zero| is time to
% construct zero constraint matrix; all times in seconds.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-22
% * Changed:    2017-11-06
%
%%

assert(size(xa, 2) == size(xb, 2), 'xa and xb must have same number of columns.');
assert(size([xa; xb],1) == size(z,1), '[xa; xb] and y must have same number of rows.');

% START measure time full computation
time.all = cputime;

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

% continuity constraint
if ~exist('x0', 'var') || ~isnumeric(x0)
    varargin = [{x0} varargin];
    x0 = NaN;
end

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

% problem structure
problem.solver = 'lsqlin';
% use active-set algorithm (depricated)
problem.options = optimoptions(problem.solver, 'Algorithm', 'interior-point');


%% Continuity constraint
% Aeq1*q1 - Aeq1*q2 = 0

% START measure time construction Aeq
time.eq = cputime;

if isnan(x0)
    % no continuity constraint
    Aeq = [];
    beq = [];
elseif m == 1
    Aeq1 = double(p(x0)');
    Aeq = [Aeq1 -Aeq1];
    beq = 0;
elseif m >= 2
    one = num2cell(ones(1,m-1));
    [~, ~, rtilde] = monomials(n, m-1);
    Aeq1 = zeros(rtilde,r);
    j = 0;
    for N=0:n
        [pN, ~, rN] = monomials({N}, m);
        pNx0 = double(pN(x0,one{:})');
        Aeq1(1:rN,j+(1:rN)) = diag(pNx0); %(rN:-1:1));
        j = j + rN;
    end
    Aeq = [Aeq1 -Aeq1];
    beq = zeros(rtilde,1);
else
    error('Continuity constraint for more than 2 variables is not supported yet.');
end

% STOP measure time construction Aeq
time.eq = cputime - time.eq;


%% Zero constraint
% Azero*q = 0

% START measure time construction Azero
time.zero = cputime;

if isempty(y0) || all(isnan(y0))
    % no constraint
    Azero = [];
    bzero = [];
else
    Azero1 = eye(r);
    if length(y0) < m
        y0 = [ones(1,m-length(y0)) y0];
    end
    Y = num2cell(y0);
    pY = double(p(Y{:}));
    Azero1(pY==0,:) = [];
    r0 = size(Azero1,1);
    Azero = [Azero1 zeros(r0,r); zeros(r0,r) Azero1];
    bzero = zeros(2*r0,1);
end

% STOP measure time construction Azero
time.zero = cputime - time.zero;


%% least squares objective
% find q minimizing the L2-norm
% ||C*q-d||^2

% START measure time construction C, d
time.obj = cputime;

problem.C = zeros(ka+kb, 2*r);
problem.d = z;
for j = 1:ka
    Xj = num2cell(xa(j,:));
    problem.C(j,1:r) = double(p(Xj{:})');
end
for j = 1:kb
    Xj = num2cell(xb(j,:));
    problem.C(ka+j, r+(1:r)) = double(p(Xj{:})');
end

% remove NaN rows from C, d
In = isnan(z);
problem.C(In,:) = [];
problem.d(In)   = [];

% STOP measure time construction C, d
time.obj = cputime - time.obj;


%% Linear least square problem
% solve LSQ min||C*q - d|| for q
% where Aineq*q <= bineq
problem.Aineq = ones(1,2*r);
problem.bineq = 1e4;
% and Aeq*q == beq
problem.Aeq = [Aeq; Azero];
problem.beq = [beq; bzero];


% START measure time solving LSQ
time.lsq = cputime;

[q, resnorm] = lsqlin(problem);

% STOP measure time solving LSQ
time.lsq = cputime - time.lsq;


%% Return fitobject & GoF
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

% RMSE is square root of residual norm
gof.rmse = sqrt(resnorm);


% STOP measure time full computation
time.all = cputime - time.all;

end
