function varargout = tolatex(obj, data_file, var, varargin)
%TOLATEX Writes latex function representation of piece-wise fit.
%
%% Usage and description
%
%   fileID = tolatex(obj, data_file, var, name)
%
% writes to file |data_file|, where
%
% * |var|   is cell array of variable name strings, default 'x';
% * |name|  is function name, default 'f';
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-01-30
% * Changed:    2018-01-30
%
%%

% number of cases
m = size(obj(1).coeffs,2);

% determine recursion
if ~isstruct(data_file)
% no case selected
j = -1;

% determine free variable
if ~exist('var', 'var') || isempty(var)
    p.var = obj(1).var;
elseif ~iscell(var)
    p.var = {var};
else
    p.var = var;
end

% determine function name
for i=1:length(varargin)
    arg = varargin{i};
    if ~isfield(p, 'dict') && iscell(arg),          p.dict = arg;
    elseif ~isfield(p, 'name') && ~isempty(arg),    p.name = arg;
    end
end
if ~isfield(p, 'name'), p.name = obj(1).name;   end
if ~isfield(p, 'dict'), p.dict = {};            end


% open file for writing
if ~isempty(data_file)
    p.file = fopen(data_file, 'w', 'n', 'UTF-8');
else
    % write to stdout
    p.file = 1;
end


p.vfmt = '\\num{%#.4e}';
p.lfmt = '\\SI{%#.4e}{\\degree}';
p.lcnv = @rad2deg;

fprintf(p.file, ...
        '%% THIS FILE HAS BEEN WRITTEN BY pwfitobject/tolatex.m %%\n\n');
    
    
if ~isempty(obj(1).xi)
    fprintf(p.file, ['\\providecommand\\%slim{' p.lfmt '}\n\n'], p.var{1}, p.lcnv(obj(1).xi));
end

tolatex(obj, p, -1);

if p.file > 2 % other than stdout/stderr
    fclose(p.file);
end

varargout = {p.file};
    
else
% first input is p-struct, second selected case
p = data_file;
j = var;

if length(obj) > 1
    for o = obj(:)'
        fprintf(p.file, '%%%% %s(%s)\n', o.name, parameter(o.var));
        
        tolatex(o, p, -1);
        
        fprintf(p.file, '\n');
    end    
elseif m > 1 && j < 0
    for j=1:m
        tolatex(obj, p, j);
    end
else
    fprintf(p.file, '\\newcommand\\%s', replace(obj.name, {'.' '_'}, ''));
    if j < 0
        j = 1;
    else
        fprintf(p.file, pad('', j, 'i'));
    end
    
    if ~isempty(obj.var)
        var = obj.var;
        
        if ~isempty(p.dict)
            [~, I] = ismember(var, p.var);
            var = p.dict(I);
        end
    else
        var = p.var;
    end
    
    tex = totex(obj, var, p.vfmt, p.lfmt, p.lcnv, [], {'^'}, ' ', [], j);
    fprintf(p.file, '{%s}\n', tex);
    
    if j < m
    end
end

end

end

function par = parameter(var)
    if isempty(var), var = {}; end
    
    par = [sprintf('%s,', var{1:end-1}) var{end}];
end