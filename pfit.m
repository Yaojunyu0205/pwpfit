function [fitobject,gof,time] = pfit (x, z, n, varargin)
%PFIT Fits multi-dimensional, polynomial function to data.
%
% Finds polynomial function in x1,...,xm of degrees n
%
%   f(x) = bn0 x1^n + ... + b0n xm^n + ... + b10 x1 + ... + b01 xm + b0;
%
% minimizing
%
%   sum[j=1:k] |f(x1(j),...,xm(j)) - z(j)|^2,
%
% where k is the length of x and z.
%
%% Usage and description
%
%   fitobject = pfit(x, z, n)
%   [...] = pfit(..., [y0 | NaN], [w | 1], [pwfoargs...])
%
% Returns fit of x against z, where x and z are column vectors with
% size(x) = size(z).
% Weights w must be a scalar (no weighting) or vector with 
% size(w) == size(z).
%
% If the optional parameter y0 (and |y0 != NaN|) is given, the returned fit
% is zero in the parameter xj if and only if the j-th component of y0 is
% zero; i.e.
%
%   f(...,xj=0,...) = 0
%
% for all x1,...,x[j-1],x[j+1],...,xm in R^(m-1).
%
% The optional argument(s) |pwfoargs| are applied to |fitobject|.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-23
% * Changed:    2018-10-10
%
%%

time.type = 'pfit';

% START measure time full computation
time.all = cputime;

% number of columns in data
m = size(x, 2);

% column of monomials to degrees n1,...,nm
% p(x) = [1,...,x^n]^T
% where the length of p is r.
[p, X, r] = monomials(n, m);

% length of data
% k = #x = #y = #z
k = length(z);

if ~isempty(varargin) && isnumeric(varargin{1})
    y0 = varargin{1};
    varargin(1) = [];
else
    y0 = NaN;
end

% weights
if ~isempty(varargin) && isnumeric(varargin{1})
    w = varargin{1};
    W = diag(w);
    varargin(1) = [];
else
    W = 1;
end

% assert number of data rows
assert(size(x,1) == size(z,1), 'x and z must have same number of rows.');

% assert number of weight elements
assert(isscalar(W) || length(W) == size(z,1), 'W and z must have same number of elements');

%% Reduction to least-square optimization
%
% As f is polynomial of degree n, i.e.
%
%   f = q0 + q10 x1 + ... + q01 xm + ... + qn0 x1^n + ... + q0n xm^n + ,
%
% with q = [q0 q10 ... q01 ... qn0 ... q0n]^T, the objective can be written
% as least-square problem in q:
%
%   find q minimizing || W*C*q - z ||^2,
%
% where ||.|| is the L2-norm,
%
%
%       | w1 |  0 |
%   W = |    \    |,
%       | 0  | wk |
%
% and
%
%       | 1 x1,1 ... xm,1 ... x1,1^n ... xm,1^n  |
%   C = | :   :   \    :   \    :     \    :     |.
%       | 1 x1,k ... xm,k ... x1,k^n ... xm,k^n  |
%
%
% If set, C is subject to the zero constraint f(Y0) == 0, 
% which for m = 1 is equivalent to the matrix equality
%
%       [1 y0 ... y0^n]*q = 0
%
% For m = 2, we have the surface equality constraint
%
%   f(x, y0) = 0
%
% for all y in R equivalent to the matrix equality
%
%       [1, ..., x^n]*A'*q = 0
%   <=>
%       A'*q = [0, ..., 0]^T
%
% with
%
%        |  1  |  ...  | 0 ... y0^(n-1) | 0 ...    0     y0^n |
%        |-----|       |     /          |       y0^(n-1)      |
%   A' = |-------------| 1              |     /               |.
%        |------------------------------| 1                 0 |
%        |----------------------------------------------------|
%

% problem structure
problem.solver = 'lsqlin';
% use active-set algorithm (depricated)
problem.options = optimoptions(problem.solver, 'Algorithm', 'active-set');


%% Zero constraint
% Aeq*q = 0

% START measure time construction Azero
time.zero = cputime;

if isempty(y0) || all(isnan(y0))
    % no zero constraint
    Azero = [];
    bzero = [];
else
    Azero = eye(r);
    if length(y0) < m
        y0 = [ones(1,m-length(y0)) y0];
    end
    Y = num2cell(y0);
    pY = double(p(Y{:}));
    Azero(pY==0,:) = [];
    r0 = size(Azero,1);
    bzero = zeros(r0,1);
end

% STOP measure time construction Azero
time.zero = cputime - time.zero;


%% least squares objective
% find q minimizing the L2-norm
% ||C*q-d||^2

% START measure time construction C, d
time.obj = cputime;

C = zeros(k, r);
for j = 1:k
    Xj = num2cell(x(j,:));
    C(j,:) = double(p(Xj{:})');
end
% ||W*(Cq - d)|| = ||W*Cq - W*d||
% for W positive diagonal
problem.C = W*C;
problem.d = diag(W).*z;

% remove NaN rows from C, d
In = isnan(z);
problem.C(In,:) = [];
problem.d(In)   = [];

% STOP measure time construction C, d
time.obj = cputime - time.obj;


%% Linear least square problem
% solve LSQ min||C*q - d|| for q
% where Aineq*q <= bineq
problem.Aineq = ones(1,r);
problem.bineq = 1e4;
% and Aeq*q == beq
problem.Aeq = Azero;
problem.beq = bzero;


% START measure time solving LSQ
time.lsq = cputime;

[q, resnorm] = lsqlin(problem);

% STOP measure time solving LSQ
time.lsq = cputime - time.lsq;


%% Return fitobject & GoF
% function
P = formula(p);
F = q'*P;
f = symfun(F, X);

fitobject = pwfitobject(['poly' sprintf('%g', n+zeros(1,m))], f, [], q, n, varargin{:});

% RMSE is square root of residual norm
gof.rmse = sqrt(resnorm);


% STOP measure time full computation
time.all = cputime - time.all;

end
