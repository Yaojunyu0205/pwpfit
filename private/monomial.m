function [P, l] = monomial(P, l, X, n, X0)
%MONOMIAL Creates monomials in variables X of degree n.

    m = length(X);
    
    if ~exist('X0', 'var')
        X0 = 1;
    end

    switch m
        case 1, P(l) = X0*X^n; l = l+1;
        otherwise
            for j = 0:n
                [P, l] = monomial(P, l, X(2:end), j, X0*X(1)^(n-j));
            end
    end
end