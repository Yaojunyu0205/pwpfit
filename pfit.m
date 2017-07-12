function [fitobject] = pfit (x, z, n, varargin)
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
% If the zero constraint is set, the solution f fulfills
%
%   f(x1,Y0) = 0
%
% for all x1.
%
%% Usage and description
%
%   fitobject = pwpfit([x1,...,xm], z, [n1,...,nm], zero)
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-23
% * Changed:    2017-02-23
%
%%

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


%% Reduction to least-square optimization
%
% As f is polynomial of degree n, i.e.
%
%   f = q0 + q10 x1 + ... + q01 xm + ... + qn0 x1^n + ... + q0n xm^n + ,
%
% with q = [q0 q10 ... q01 ... qn0 ... q0n]^T, the objective can be written
% as least-square problem in q:
%
%   find q minimizing || C*q - z ||^2,
%
% where ||.|| is the L2-norm and
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


% zero constraint
% Aeq*q = 0
if isempty(y0) || all(isnan(y0))
    % no zero constraint
    Azero = [];
    bzero = [];
elseif m == 1
    Azero = double(p(y0)');
    bzero = 0;
elseif m == 2
    Azero = zeros(n+1,r);
    j = 0;
    for N=0:n
        [pN, ~, rN] = monomials(N, 1);
        pNy0 = double(pN(y0)');
        Azero(1:rN,j+(1:rN)) = cdiag(pNy0(rN:-1:1));
        j = j + rN;
    end
    bzero = zeros(n+1,1);
elseif m > 2
    Azero = zeros(n+1,r);
    Y0 = num2cell(y0);
    j = 0;
    for N=0:n
        for i=0:N
            % get monomials vector in x2,...,xm of degree i
            [pNi, ~, rNi] = monomials({i}, m-1);
            pNiy0 = double(pNi(Y0{:}));
            Azero(1+N-i,j+(1:rNi)) = pNiy0';
            j = j + rNi;
        end
    end
    bzero = zeros(n+1,1);
else
    error('Zero constraint for more than 2 variables is not supported yet.');
end


% least squares objective
% find q minimizing the L2-norm
% ||C*q-d||^2
C = zeros(k, r);
d = z;
for j = 1:k
    Xj = num2cell(x(j,:));
    C(j,:) = double(p(Xj{:})');
end


% inequality condition
% A*q <= b
A = ones(1,r);
b = 1e4;

% solve LSQ for q
q = lsqlin(C, d, A, b, Azero, bzero);

% function
P = formula(p);
F = q'*P;
f = symfun(F, X);

fitobject = pwfitobject(['poly' sprintf('%g', n+zeros(1,m))], f, [], q, n, varargin{:});


end
