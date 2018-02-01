% Piece-wise identification of the GTM aerodynamic coefficients.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2017-11-01
% * Changed:    2017-11-20
%
%% Variables, constants, and their units
%
% * |alpha|  :  angle of attack,                                    rad
% * |beta|   :  side-slip angle,                                    rad
% * |zeta|   :  rudder deflection,                                  rad
% * |eta|    :  elevator deflection,                                rad
% * |xi|     :  aileron deflection,                                 rad
% * |Cl|     :  aerodynamic coefficient moment body x-axis,         -
% * |Cm|     :  aerodynamic coefficient moment body y-axis,         -
% * |Cn|     :  aerodynamic coefficient moment body z-axis,         - 
% * |Cx|     :  aerodynamic coefficient force body x-axis,          -
% * |Cy|     :  aerodynamic coefficient force body y-axis,          -
% * |Cz|     :  aerodynamic coefficient force body z-axis,          -
% * |phat|   :  normalized roll rate,                               rad
% * |qhat|   :  normalized pitch rate,                              rad
% * |rhat|   :  normalized yaw rate,                                rad
%%


%% GTM aerodynamic coefficients
% load 6-DOF aerodynamics data 
% with inputs |p|, aerodynamic coefficients |Aero|
load('aerodata.mat');


% pre-stall data boundary
istar = 14;
alphastar = p.alpha(istar);


%% Find stall angle of attack
% fit Cx with respect to alpha
% to get |alpha0| such that
%
%   CxPre(alpha0) == CxPost(alpha0).
[Cx1, alpha0, ~, time.x1] = pwpfit(p.alpha(1:istar), p.alpha(istar+1:end), ...
                     Aero.Cx(:,p.beta==0,p.xi==0,p.eta==0,p.zeta==0), ...
                     3, NaN);

%% Plot resulting fit
figure(1)
plot(Cx1, deg2rad([-10 90]))
hold on
% plot GTM data for comparison
plot(p.alpha, Aero.Cx(:,p.beta==0,p.xi==0,p.eta==0,p.zeta==0), '+')
hold off


%% Prepare data for fitting
% aerodynamic data is given as M-dimensional matrizes
% with dimensions [alpha]x[beta]x[xi]x[eta]x[zeta].

% prepare data to obtain
%   ( alpha(i),beta(i),xi(i),eta(i),zeta(i),
%       CX(i),CY(i),CZ(i),CL(i),CM(i),CN(i) )
% for 1 <= i <= k.
[alpha, beta, xi, eta, zeta, CX] = prepareHyperSurfaceData(p.alpha, p.beta, p.xi, p.eta, p.zeta, permute(Aero.Cx,5:-1:1));
[  ~  ,  ~  , ~ ,  ~ ,  ~  , CY] = prepareHyperSurfaceData(p.alpha, p.beta, p.xi, p.eta, p.zeta, permute(Aero.Cy,5:-1:1));
[  ~  ,  ~  , ~ ,  ~ ,  ~  , CZ] = prepareHyperSurfaceData(p.alpha, p.beta, p.xi, p.eta, p.zeta, permute(Aero.Cz,5:-1:1));
[  ~  ,  ~  , ~ ,  ~ ,  ~  , CL] = prepareHyperSurfaceData(p.alpha, p.beta, p.xi, p.eta, p.zeta, permute(Aero.Cl,5:-1:1));
[  ~  ,  ~  , ~ ,  ~ ,  ~  , CM] = prepareHyperSurfaceData(p.alpha, p.beta, p.xi, p.eta, p.zeta, permute(Aero.Cm,5:-1:1));
[  ~  ,  ~  , ~ ,  ~ ,  ~  , CN] = prepareHyperSurfaceData(p.alpha, p.beta, p.xi, p.eta, p.zeta, permute(Aero.Cn,5:-1:1));
            
% data input vector
X = [alpha beta xi eta zeta];


%% Fit longitudinal coefficients
% as 5-D polynomials in X
% with constraint of continuity in |alpha0|
[Cx, ~, gofx, time.x] = pwpfit(X(alpha<=alphastar,:), X(alpha>alphastar,:), CX, 3, alpha0, 'Cx');
[Cz, ~, gofz, time.z] = pwpfit(X(alpha<=alphastar,:), X(alpha>alphastar,:), CZ, 3, alpha0, 'Cz');
[Cm, ~, gofm, time.m] = pwpfit(X(alpha<=alphastar,:), X(alpha>alphastar,:), CM, 3, alpha0, 'Cm');


%% Fit lateral coefficients
% as 5-D polynomials in X
% with constraint of continuity in |alpha0| and
% with zero constraint for beta=0, xi=0, zeta=0
[Cy, ~, gofy, time.y] = pwpfit(X(alpha<=alphastar,:), X(alpha>alphastar,:), CY, 3, alpha0, [1 0 0 1 0], 'Cy');
[Cl, ~, gofl, time.l] = pwpfit(X(alpha<=alphastar,:), X(alpha>alphastar,:), CL, 3, alpha0, [1 0 0 1 0], 'Cl');
[Cn, ~, gofn, time.n] = pwpfit(X(alpha<=alphastar,:), X(alpha>alphastar,:), CN, 3, alpha0, [1 0 0 1 0], 'Cn');


%% Export coefficients
% write coefficients to file
tomatlab([Cx Cy Cz Cl Cm Cn], 'coefficients.m', ...
         {'alpha' 'beta' 'xi' 'eta' 'zeta'})
