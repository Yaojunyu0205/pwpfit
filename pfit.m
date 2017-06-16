function [fitobject, x0] = pfit (x, z, n)
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
%   fitobject = pwpfit([x1,...,xm], z, [n1,...,nm])
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-23
% * Changed:    2017-02-23
%
%%

% column of monomials to degrees n1,...,nm
% p(x) = [1,...,x^n]^T
% where the length of p is r.
[p, X, r] = monomials(n, length(x(1,:)));

% length of data
% k = #x = #y = #z
k = length(z);


%% Reduction to least-square optimization
%
% As f is polynomial of degree n, i.e.
%
%   f = qn0 x1^n + ... + q0n xm^n + ... + q10 x1 + ... + q01 xm + q0,
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
q = lsqlin(C, d, A, b);

% function
P = formula(p);
F = q'*P;
f = symfun(F, X);

fitobject = pwfitobject(sprintf('poly%g%g', n, n), f, [], q, n);


end
