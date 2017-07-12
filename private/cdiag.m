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
