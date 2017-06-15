function tex = totex(obj, var, vfmt, lfmt, lcnv, order, efmt, mfmt)
%TOTEX Returns tex string representation of piece-wise fit.
%
%% Usage and description
%
%   tex = totex(obj, var, vfmt, lfmt, lcnv, order, efmt, mfmt)
%
% where
%
% * |var|   is cell array of variable name strings, default 'x';
% * |vfmt|  is format spec for coefficients values, default '%#.4e', or
%               cell array of left- and right-hand side delimiter, 
%               e.g. {'(' ')'} equals to '(%#.4e)';
% * |lfmt|  is format spec for limit values, default is |vfmt|;
% * |lcnv|  is conversion function for limit values, default is |id|;
% * |order| is order of monomials, either 'asc' or 'desc', default is
%               ascending by degree;
% * |efmt|  is format spec for exponentials, default is '^{%u}', or cell
%               array {[|exp|] [|efl| |efr|]} equals to [exp efl '%u' efr]
%               with defaults are exp='^', efl='', efr='';
% * |mfmt|  is format spec for multiplication sign, default is ' '.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-27
% * Changed:    2017-06-14
%
%%

% number of cases
m = size(obj.coeffs,2);

% determine recursion
if nargin > 1 && isstruct(var)
    % first input is p-struct, second selected case
    p = var;
    j = vfmt;
else
% case selected
j = -1;

% determine free variable
if nargin < 2 || isempty(var)
    p.var = 'x';
else
    p.var = var;
end

% determine value-to-tex function
if nargin < 3 || isempty(vfmt)
    p.vfmt = '%#.4e';
elseif iscell(vfmt)
    p.vfmt = [vfmt{1} '%#.4e' vfmt{2}];
else
    p.vfmt = vfmt;
end

% determine limit-to-tex function
if nargin < 4 || isempty(lfmt)
    p.lfmt = p.vfmt;
elseif iscell(lfmt)
    p.lfmt = [lfmt{1} '%#.4e' lfmt{2}];
else
    p.lfmt = lfmt;
end

% determine limit conversion function
if nargin < 5 || isempty(lcnv)
    p.lcnv = @(x) x;    %id
else
    p.lcnv = lcnv;
end

% determine order of terms
if nargin < 6 || isempty(order) || strcmp(order, 'asc')
    p.order   = 'asc';
    p.degrees = 0:obj.degree;
else
    p.order   = 'dsc';
    p.degrees = obj.degree:-1:0;
end

% determine exponential-to-tex function
if nargin < 7 || isempty(efmt)
    p.efmt = '^{%u}';
elseif iscell(efmt)
    if length(efmt) == 1,   efmt = [efmt {'' ''}];
    elseif length(efmt) == 2, efmt = [{'^'} efmt];
    end
    
    p.efmt = [efmt{1} efmt{2} '%u' efmt{3}];
else
    p.efmt = efmt;
end

% determine multiplication symbol
if nargin < 8 || isempty(mfmt)
    p.mfmt = ' ';
else
    p.mfmt = mfmt;
end

end

tex = '';
if m > 1 && j < 0
    tex = sprintf('%s\\left\\{\\begin{array}{c l} \n', tex);
    for j=1:m
        tex = sprintf('%s%s', tex, totex(obj, p, j));
        if j < m, tex = sprintf('%s \\\\\n', tex); end
    end
    tex = sprintf('%s\t& \\text{else}\n\\end{array}\\right.', tex);
else
    if j < 0
        j = 1;
    end
    l = 1;
    for i=p.degrees
        if ~iscell(p.var)
            tex = printterm(tex, p, obj.coeffs(l,j), i);
            l = l + 1;
        else
            for k=0:i
                tex = printterm(tex, p, obj.coeffs(l,j), {i k});
                l = l + 1;
            end
        end
    end
    if j < m 
        tex = sprintf(['%s\t& \\text{if $%s \\leq ' p.lfmt '$}'], tex, p.var, p.lcnv(obj.xi(j)));
    end
end


end

function tex = printterm(tex, p, coeff, i)
    if coeff >= 0 && ~strcmp(tex, '')
        tex = sprintf('%s+ ', tex);
    elseif coeff < 0
        tex = sprintf('%s- ', tex);
    end
    tex = sprintf(['%s' p.vfmt '%s '], tex, abs(coeff), monomial(p, i));
end
    
function tex = monomial(p, i)
    if ~iscell(p.var)
%         if iscell(i), i = i{1}; end
        switch i
            case 0, tex = '';
            case 1, tex = sprintf([p.mfmt '%s'], p.var);
            otherwise
                tex = sprintf([p.mfmt '%s' p.efmt ''], p.var, i);
        end
    else
        tex = [monomial(getvar(p,1), i{1}-i{2}) monomial(getvar(p,2), i{2})];
    end
end

function pout = getvar(pin, k)
    %GETVAR     Returns p-struct for k-th variable.
    
    pout = pin;
    pout.var = pin.var{k};
end
