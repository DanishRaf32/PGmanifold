function [yROM,xxROM,Quadratic_ROM_time,deltar,omegar]=QuadROM_solve(f,Jack,tspan,x0,solver,V,n_oc,Vbar,W)
xdotr= @(t,x) V'*f(V*x + Vbar*get_x_sq(x')'); %ode
x0r=V'*x0;
optionsROM = odeset('RelTol',1e-6,'AbsTol',1e-9);
optionsNL.solver = 'fsolve';   %'fsolve','NewtonRaphson'
optionsNL.funJacobian=Jack;
optionsNL.RelTol = 1e-6;
optionsNL.AbsTol = 1e-9;
tic
switch(solver)
    case 'ode15s'
       [tFOM,xROM]=ode15s(xdotr,tspan,x0r,optionsROM);
    case 'ode45'
       [tFOM,xROM]=ode45(xdotr,tspan,x0r,optionsROM);
    case 'Euler'
        [tFOM,xROM]=implicitEuler(xdotr,tspan,x0r,optionsNL);
end
Quadratic_ROM_time=toc
xxROM = V*xROM' + Vbar*W; % check the dimension = (n x ts)
%xxROM=(V*xROM')';
deltar=xxROM(1:n_oc,:);
omegar=xxROM(n_oc+1:end,:);
avg_delr=sum(deltar,1)/n_oc; % 1st state-variable
avg_omegar=sum(omegar,1)/n_oc; %2st state-variable 
yROM(:,1)=avg_delr;
yROM(:,2)=avg_omegar;
end