function varargout = plot(obj, varargin)
    %PLOT Plots the piece-wise defined function.
    %
    % See also PLOT
    %%
    
    X = symvar(obj.f);
	m = length(X);
    
    switch m
        case 2, varargout{:} = fsurf(obj.f, varargin{:});
            
        case 1, varargout{:} = fplot(obj.f, varargin{:});
        
        otherwise
            b = ceil(sqrt(m-1));
            a = ceil((m-1)/b);
            
            ij = 0;
            varargout = cell(1,m-1);
            
            for i=1:a
                for j=1:b
                    ij = ij + 1;
                    if ij > m-1, break; end
                    
                    %else
                    subplot(a, b, ij);
                    Xab = num2cell(X);
                    Xab((0:m-1)~=0 & (0:m-1)~=ij) = {0};
                    varargout{ij} = fsurf(obj.f(Xab{:}), varargin{:});
                end
            end
    end
                

    if nargout < 1
        varargout = {};
    end
end