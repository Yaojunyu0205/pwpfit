classdef pwfitobject
%PWFITOBJECT Result of piece-wise fitting.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-02-23
% * Changed:    2017-02-23
%
%% Variables
%
% * |f|      :  fitted, piece-wise defined function
% * |fi|     :  fitted sub-functions
% * |m|      :  number of sub-functions, i.e. m = #fi = #xi
% * |x|      :  variable(s) of f
% * |xi|     :  upper boundary for fi
%
%%

    properties
        type;
        f;
        
        fi;
        xi;
        
        degree;
        coeffs;
    end
    
    methods(Static, Access=protected)
        function f = pwfunction(fi, xi)
           %PWFUNCTION Returns a piece-wise defined symbolic function.
           %%
           m = length(fi);
           assert(m == length(xi)+1);
           
           X = symvar(fi{1});
           
           varargs = cell(1,2*m-1);
           for i = 1:m-1
               varargs{i}   = (X(1) <= xi(i));
               varargs{i+1} = fi{i};
           end
           varargs{end} = fi{end};
           
           f(X) = piecewise(varargs{:});
        end
    end
    
    methods
        tex = totex(obj, var, vfmt, lfmt, lcnv, order, efmt, mfmt);
        
        function obj = pwfitobject(type, fi, xi, coeffs, degree)
            %PWFITOBJECT Creates a new pwfitobject.
            %
            %% Usage and description
            % 
            %   obj = pwfitobject(type, fi, xi, coeffs)
            %
            % Where the number of sub-functions |fi| must equal the number 
            % of columns in |coeffs|, and must exceed the number of limits
            % |xi| by one, |#fi = #xi + 1|.
            %%
            if ~iscell(fi), fi = {fi}; end

            obj.type = type;
            
            obj.fi = fi;
            obj.xi = xi;
            obj.coeffs = coeffs;
            obj.degree = degree;
            
            obj.f = pwfitobject.pwfunction(fi, xi);
        end
        
        function varargout = plot(obj, varargin)
            %PLOT Plots the piece-wise defined function.
            %
            % See also PLOT
            %%
            if length(symvar(obj.f)) > 1
                varargout{:} = fsurf(obj.f, varargin{:});
            else
                varargout{:} = fplot(obj.f, varargin{:});
            end
            
            if nargout < 1
                varargout = {};
            end
        end
        
        function f = plus(obj1, obj2)
            %PLUS Overloaded binary operator '+'.
            %   Returns the piece-wise defined function f = f1 + f2, where
            %   f1, f2 are the piece-wise fitted functions of obj1, obj2.
            m1 = length(symvar(obj1.f));
            m2 = length(symvar(obj2.f));
            x = sym('X', [max(m1,m2) 1]);
            X = num2cell(x);
            f(x) = obj1.f(X{1:m1}) + obj2.f(X{1:m2});
        end
        
        function varargout = subsref(obj, s)
            %SUBSREF Overloaded subsref function.
            %
            % See also SUBSREF
            %%
            switch s(1).type
                case '()'
                    varargout = {double(obj.f(s.subs{:}))};
                otherwise
                    varargout{:} = builtin('subsref',obj,s);
            end
        end
    end
end

