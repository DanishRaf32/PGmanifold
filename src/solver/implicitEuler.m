function [tsim,y,nNLSE] = implicitEuler(odefun,tsim,y0,varargin)
% IMPLICITEULER - Implicit (backward) Euler scheme for simulation
%
% Syntax:
%		[tsim,y] = IMPLICITEULER(odefun,tsim,y0)
%		[tsim,y] = IMPLICITEULER(odefun,tsim,y0,Opts)
%
% Description:
%       This function solves the system odefun over the time interval tsim 
%       and with the initial condition y0, using an implicit (backward) Euler
%       method. 
%
%       For the solution of the resulting nonlinear systems of equations (NLSE),
%       either MATLAB's built-in function fsolve or the function NewtonRaphson
%       can be selected.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -odefun:    right hand-side of system ODE, yDot = odefun(t,y)
%       -tsim:      simulation time vector s.t. tsim = t0:h:tend
%       -y0:        initial condition for y
%       *Optional Input Arguments:*
%       -Opts:         a structure containing following options
%           -.funJacobian:
%           -.solver:
%           -.MaxIter:
%           -.Display:
%           -.TolX:
%           -.TolFun:
%           -.Algorithm:
%           -.SpecifyObjectiveGradient:
%
% Output Arguments:
%       -tsim:  simulation time vector s.t. tsim = t0:h:tend
%       -y:     solution matrix, rows are solution states yi
%
% Examples:
%
% See Also: 
%       NewtonRaphson, fsolve
%
% References:
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch pod-deim">pod-deim</a>, a   
% Simulation-Based Model Order Reduction Toolbox based on POD-DEIM developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch pod-deim">pod-deim</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Maria Cruz Varona, Julian Suk, Nico Schneucker 
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  18 Jul 2018
% Copyright (c) 2016-2018 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% Parse options
if nargin >= 4
    Opts = varargin{1};  
end
% Options for implicitEuler
Def.solver = 'NewtonRaphson'; % set solver to 'NewtonRaphson' or 'fsolve'
Def.funJacobian = []; %[] means numerical Jacobian by finite-difference

% Options for 'fsolve' AND 'NewtonRaphson'
Def.MaxIter = 400;
Def.Display = 'none'; % 'none' or 'iter'

% Options for 'fsolve'
Def.TolX = 1e-6;
Def.TolFun = 1e-6;
Def.Algorithm = 'trust-region-dogleg'; % 'trust-region-dogleg', 'trust-region-reflective', 'levenberg-marquardt'
Def.SpecifyObjectiveGradient = true;

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%%
% number of NLSE
nNLSE = length(tsim);

% Dimension of ode
n = length(y0);

% Uniform step size h
h = tsim(2)-tsim(1);

% Initial time step and y
tnew = tsim(1) + h;
ycurr = y0;

% Preallocation of y
y = zeros(length(tsim),length(y0));
y(1,:) = y0;

% Numerical or analytical Jacobian:
if isempty(Opts.funJacobian)
    FJacTimeDep = @(t, x) -speye(n) + h*NumJacobian(@(a) odefun(t, a),x,Opts);
else
    FJacTimeDep = @(t, x) -speye(n) + h*Opts.funJacobian(t,x);
end
    
% Implicit Euler method 
for i = 1:length(tsim)-1
    % if we want to substitute Euler with another quadrature formula, only the following line would have to change
    sys = @(ynew) ycurr - ynew + h*odefun(tnew,ynew);
    FJac = @(x) FJacTimeDep(tnew, x);
    
    if (strcmp(Opts.solver,'fsolve'))
        [ynew] = fsolve(inputfsolve(sys,FJac),ycurr,Opts);
    else
        [ynew] = NewtonRaphson(sys,ycurr,FJac,Opts);
    end
    y(i+1,:) = ynew;
    tnew = tnew + h;
    ycurr = ynew;
end
end

function [sys,FJac] = inputfsolve(sys,FJac)
end