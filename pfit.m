function [fitobject, x0] = pfit (x, z, n, varargin)
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
%   f(x0,x2,...,xm) = 0
%
% for all [x2...xm].
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
    x0 = varargin{1};
    varargin(1) = [];
else
    x0 = NaN;
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


% zero equality constraint
% Aeq*q = 0
if isempty(x0) || isnan(x0)
    % no equality constraint
    Aeq = [];
    beq = [];
    x0 = NaN;
elseif m == 1
    Aeq = double(p(x0)');
    beq = 0;
elseif m == 2
    Aeq = zeros(n+1,r);
    j = 0;
    for N=0:n
        [pN, ~, rN] = monomials(N, 1);
        pNx0 = double(pN(x0)');
        Aeq(1:rN,j+(1:rN)) = cdiag(pNx0(rN:-1:1));
        j = j + rN;
    end
    beq = zeros(n+1,1);
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
q = lsqlin(C, d, A, b, Aeq, beq);

% function
P = formula(p);
F = q'*P;
f = symfun(F, X);

fitobject = pwfitobject(['poly' sprintf('%g', n+zeros(1,m))], f, [], q, n, varargin{:});


end

function D = cdiag(v)
%CDIAG  Returns a square matrix which counter-diagonal is v.
%
% A counter-diagonal matrix is 
%
%       | 0   n |
%   D = |   /   |
%       | 1   0 |

    A = diag(v);
    D = A(:,end:-1:1);
    
end
