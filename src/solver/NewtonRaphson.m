function [xend] = NewtonRaphson(f,x0,FJac,Opts)
% NEWTONRAPHSON - Newton-Raphson method for finding the solution for 0 = f(x)
%
% Syntax:
%       xend = NEWTONRAPHSON(f,x0,Jac)
%       xend = NEWTONRAPHSON(f,x0,Jac,Opts)
% 
% Description:
%       This function performs the well-known Newton-Raphson algorithm for
%       finding the root x of the function f(x), starting from the initial 
%       guess x0. In other words, this function tries to iteratively solve 
%       the nonlinear system of equations 0 = f(x).
%
%       The algorithm is repeated until a sufficiently accurate value is 
%       reached, i.e. until norm(f(x)) <= tol, or until the maximum number
%       of iterations is achieved.
%
%       Options contain settings for the relative and absolute tolerance, 
%       the maximum number of iterations and desired displaying information.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -f:    right hand-side of system ODE, yDot = odefun(t,y)
%       -x0:      simulation time vector s.t. tsim = t0:h:tend
%       -FJac:        initial condition for y
%       *Optional Input Arguments:*
%       -Opts:         a structure containing following options
%           -.RelTol:
%           -.AbsTol:
%           -.MaxIter:
%           -.Display:
%
% Output Arguments:
%       -xend:  simulation time vector s.t. tsim = t0:h:tend
%
% Examples:
%
% See Also: 
%       fsolve, implicitEuler, nlmm
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

%% Parse options for Newton-Raphson tolerance
Def.RelTol = 1e-6; % relative tolerance for residual error: norm(f(xcurr))
Def.AbsTol = 1e-9; % absolute tolerance for residual error: norm(f(xcurr))
Def.MaxIter = 400; % maximal number of iterations
Def.Display = 'none';
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end
%% Newton-Raphson method
xcurr = x0;
iter = 0;
tol = Opts.AbsTol;
fcurr = f(xcurr);
while norm(fcurr) > tol
    fcurr = f(xcurr);
    dxcurr = FJac(xcurr)\fcurr;
    
    xcurr = xcurr - dxcurr;
    iter = iter + 1;

    if(strcmp(Opts.Display,'iter'))
    disp([num2str(iter), ' iterations', '  norm(f(xcurr)): ', num2str(norm(fcurr))])
    end
    
    tol = Opts.RelTol*norm(fcurr)+Opts.AbsTol;
    if iter>=Opts.MaxIter
        fprintf(2,'maximum number of iterations reached without converging\n')
        fprintf(2,['norm(f(xcurr)): ' num2str(norm(fcurr)) '\n'])
        break
    end
end

xend = xcurr;
%if(strcmp(Opts.Display,'iter'))
 %   disp([num2str(iter) ' iterations'])
%end
end