function [p, X, r] = monomials(deg, m)
%MONOMIALS  Creates monomials in multiple variables.
%
%% Usage and description
%
%   [p, X, r] = monomials(n)
%
% Creates a column vector of monomials to degree n.
%
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

if ~exist('X', 'var')
    X = sym('X', [m 1]);
end

if ~iscell(deg)
    degrees = 0:deg;
    n = deg;
else
    degrees = cell2mat(deg);
    n = max(degrees);
end


% length of p for 0 < m < 3
r = ((n+1)^m + n+1)/2;
P = sym('P', [r 1]);

l = 1;
for i = degrees
    [P, l] = monomial(P, l, X, i);
end

r = l-1;
p = symfun(P(1:r), X);

end
