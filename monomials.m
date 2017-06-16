function [p, X, r] = monomials(n, m)
%MONOMIALS  Creates monomials in multiple variables.
%
%% Usage and description
%
%   [p, X, r] = monomials(n)
%
% Creates a column vector of monomials to degree n.
%
%   [p, X, r] = monomials(n, m)
%
% Creates a column vector of monomials in m variables to degree n. 
%
% p is the symbolic vector function 
%
%   p(X) = [1, x1, ..., xm, ..., x1^n, ..., xm^n]^T
%
% with X = [x1,...,xm] and p has length r.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-06-16
% * Changed:    2017-06-16
%
%%

if ~exist('m', 'var')
    m = 1;
end


% length of p for 0 < m < 3
r = ((n+1)^m + n+1)/2;
X = sym('X', [m 1]);
P = sym('P', [r 1]);
l = 1;
for i = 0:n
    switch m
        case 1, P(l) = X^i; l = l+1;
        case 2
            for j = 0:i
                P(l) = X(1)^(i-j)*X(2)^j;
                l = l+1;
            end
        otherwise
            error('Monomials of more than 2 variables are not supported yet.');
    end
end
p = symfun(P, X);

end