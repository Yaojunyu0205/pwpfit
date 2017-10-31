function [pi, piInv] = LMI2single(a)
%LMI2SINGLE Transforms a linear matrix inequality into a single inequality.
%
% Given phi(x) = a^T x < x0 with a != 0, there is a linear bijection
% pi(y) = x such that
%
%   phi(pi(y)) = y1
%
% for any x = (x1,...,xM), y = (y1,...,yM).
%
%% Usage and description
%
%   [pi, piInv] = LMI2single(a)
%
% Returns the bijection |pi| for phi(x) = a^T x such that phi(pi(y)) = y1;
% and its inverse |pi^(-1)|.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-10-30
% * Changed:    2017-10-31
%
%%

assert(a(1) ~= 0, 'The first component of |a| must be non-zero.');

piInv = eye(length(a));

% for k=1:length(a)
%     piInv(k,k:end) = a(1:end+1-k);
% end

piInv(1,:) = a;

pi = inv(piInv);