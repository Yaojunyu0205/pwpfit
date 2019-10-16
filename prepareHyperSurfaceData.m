function varargout = prepareHyperSurfaceData(varargin)
%PREPAREHYPERSURFACEDATA Prepares data for hyper-surface fitting.
%
% See also PREPARECURVEDATA, PREPARESURFACEDATA.
%
%% Usage and description
%
%   [X1out,...,XNout,Zout] = prepareHyperSurfaceData(X1in,...,XNin, Zin)
%
% Transforms data for fitting with the PFIT / PWPFIT function. |X1in|, ...,
% |XNin| are vectors of length k1,...,kN; |Zin| is an N-D matrix with
% dimensions kNx...xk1.
%
% Returns vectors |X1out|, ..., |XNout|, and |Zout| of equal length
% (k1*...*kN).
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-07-07
% * Changed:    2017-11-01
%
%%


assert(nargin == nargout, 'Number of inputs must equal number of outputs.');

varargout = cell(1, nargout);


switch(nargin)
    case 1
        %% 1-Dimension
        varargout{1} = varargin{1};
        
%     case 2
%         %% 2-Dimension: curve fitting
%         varargout{:} = prepareCurveData(varargin{:});
%         
%     case 3
%         %% 3-Dimension: surface fitting
%         [varargout{:}] = prepareSurfaceData(varargin{:});
        
%     case 4
    otherwise
        %% 4-Dimension: hyper-surface fitting
        % cell-vector of inputs
        Xin = varargin(1:end-1);
        % output matrix
        Zin = varargin{end};
        
        % assert dimension of Zin is kNx...xk1
        for i=1:length(Xin)
            assert(length(Xin{i})==size(Zin,length(Xin)-i+1), ...
                   'Length of input %u must equal (N-%u)-th dimension of Zin', i, i-1);
        end
        
        % re-order for mesh
        % dimensions of X1mesh,...,XNmesh must equal dimensions of Zin
        Xmesh = cell(1, length(Xin));
        [Xmesh{:}] = ndgrid(Xin{end:-1:1});
        
        % transform to columns
        data = cellfun( @(c) c(:), [Xmesh, Zin], 'UniformOutput', false );
        
        % re-re-order for output
        varargout = [data((end-1):-1:1), data(end)];
        
%     otherwise
%         error('Hyper-surfaces of dimension greater than 4 not supported yet.');
end

