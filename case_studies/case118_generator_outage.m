%Nonlinear Model Order Reduction of Power GRid Networks using Quadratic
%Manifolds
% Farhana and Danish
% Department of Electrical Engineering, IOT, University of Kashmir, India

% Simulates IEEE 118 Power System Models via MATPOWER and pg_sync_models using
% Euler. Reduction via POD and Quadratic Manifold.
% TEST CASE I: Generator Outage at 89th Bus
% Note: The generator dampings are to be changed to 10pu. 
%   References
%   [1] T. Nishikawa and A. E. Motter, Comparative analysis of existing
%       models for power-grid synchronization, New J. Phys. 17, 015012 (2015).
%
%   [2] A. E. Motter, S. A. Myers, M. Anghel, and T. Nishikawa, Spontaneous
%       synchrony in power-grid networks, Nat. Phys. 9, 191-197 (2013).
%
%   [3] Tobias, Sara Grundel, Nonlinear model reduction of dynamical power 
%       grid models using quadratization and balanced truncation (2021)
%
%   [4] Rafiq, Bazaz, Nonlinear Model Order Reduction via nonlinear moment-
%       matching with Dynamic Mode Decomposition (2020).
%% Pre-requisutes: Initialize the IEEE Power System and model configuration
clear;clc;close all
global data data2
tspan1 = 0:0.01:3;    % Pre-fault simulation time
tspan2 = 3:0.01:8;  % 1000 ms self-clearing Fault time
%pre-fault 
mpc=case118;   
mpc.ref_freq=60; 
[data,details]=EN_model(mpc); 
%fault
mpc2=case118_fault; %load increases to 700 
mpc2.ref_freq=60;
[data2,details2]=EN_model(mpc2);  
%System-Parameters
n_oc=length(data.H);   
u=@(t) 1;
solver='Euler'; % Euler, ode15s
%% Pre-Fault Conditions
disp('Simulating pre-fault conditions...')
% FOM
f1= @(x) power_func(x); % nonlinear mapping
x01=zeros(2*n_oc,1);
Jack1=@(x) Power_Jack(x);
[yFOM1,xFOM1,FOM_time1,delta1,omega1]=FOM_solve(f1,Jack1,tspan1,x01,solver,n_oc);
% POD
rdefl1=8;
[Vpod1]=basisRed(xFOM1',rdefl1);
[yROM_POD1,~,ROM_POD_time1,delta_POD1,omega_POD1]=ROM_solve(f1,Jack1,tspan1,x01,solver,Vpod1,n_oc);
% Quadratic Manifold
lambda1=1e-1; %1e-1=118, 1=300
shat1 = Vpod1'*(xFOM1');
W1 = get_x_sq(shat1')';
err1 = xFOM1'-Vpod1*shat1;
[q,~] = size(xFOM1'); 
[p,~]= size(W1);
Aplus = [W1'; sqrt(lambda1)*eye(p)];
bplus = [err1'; zeros(p,q)];
Vbar1 = (Aplus\bplus)'; %state-dependent basis
[yROM_Quad1,xROM_Quad1,QuadROM1_time,Deltar_Quad1,Omegar_Quad1]=QuadROM_solve(f1,Jack1,tspan1,x01,solver,Vpod1,n_oc,Vbar1,W1);

%%  Fault-Condition (Genearator outage at 89th bus)
disp('Simulating fault conditions...')
%FOM
f2= @(x) power_func2(x);
Jack2=@(x) Power_Jack2(x);
x02 = xFOM1(end,:)';
[yFOM2,xFOM2,FOM_time2,delta2,omega2]=FOM_solve(f2,Jack2,tspan2,x02,solver,n_oc);
%POD
rdefl2=rdefl1;
[Vpod2]=basisRed(xFOM2',rdefl2); %rdefl=30
[yROM_POD2,~,ROM_POD_time2,delta_POD2,omega_POD2]=ROM_solve(f2,Jack2,tspan2,x02,solver,Vpod2,n_oc);
% Quadratic Manifold
lambda2=1e-1; %1e-1=118, 1=300
shat2 = Vpod2'*(xFOM2');
W2 = get_x_sq(shat2')';
err2 = xFOM2'-Vpod2*shat2;
[q,~] = size(xFOM2'); 
[p,~]= size(W2);
Aplus = [W2'; sqrt(lambda2)*eye(p)];
bplus = [err2'; zeros(p,q)];
Vbar2 = (Aplus\bplus)'; %state-dependent basis
[yROM_Quad2,xROM_Quad2,QuadROM2_time,Deltar_Quad2,Omegar_Quad2]=QuadROM_solve(f2,Jack2,tspan2,x02,solver,Vpod2,n_oc,Vbar2,W2);
%% combining outputs
XFOM=[xFOM1' xFOM2'];
YFOM=[yFOM1; yFOM2];
OMEGA=[omega1;omega2];
DELTA=[delta1;delta2];
TSPAN=[tspan1 tspan2];
FOM_TIME=FOM_time1+FOM_time2;
%POD
YPodROM=[yROM_POD1; yROM_POD2];
PODROM_TIME = ROM_POD_time1+ROM_POD_time2;
Poddeltar=[delta_POD1;delta_POD2];
PodOmegar=[omega_POD1;omega_POD2];
%Quad
yQuadROM=[yROM_Quad1; yROM_Quad2];
QuadROM_time = QuadROM1_time+QuadROM2_time;
Quaddeltar=[Deltar_Quad1 Deltar_Quad2];
QuadOmegar=[Omegar_Quad1 Omegar_Quad2];
%% Errors
POD_DElerror=abs(YFOM(:,1) - YPodROM(:,1));
POD_Omgerror=abs(YFOM(:,2) - YPodROM(:,2));
Quad_DElerror=abs(YFOM(:,1) - yQuadROM(:,1));
Quad_Omgerror=abs(YFOM(:,2) - yQuadROM(:,2));
%Norms
err_PODdeltar=norm(abs(YFOM(:,1) - YPodROM(:,1)))./norm(YFOM(:,1))
err_Quaddeltar=norm(abs(YFOM(:,1) - yQuadROM(:,1)))./norm(YFOM(:,1))
err_PODomegar=norm(abs(YFOM(:,2) - YPodROM(:,2)))./norm(YFOM(:,2))
err_Quadomegar=norm(abs(YFOM(:,2) - yQuadROM(:,2)))./norm(YFOM(:,2))
%% Plotting
%Output
plotOut(TSPAN,YFOM,YPodROM,yQuadROM)
%% Error
  figure(2)
plot(122)
semilogy(TSPAN,POD_Omgerror,TSPAN,Quad_Omgerror)
legend('Linear Manifold','Quadratic Manifold')
title('Error in Omega')
grid on
% Plotting Nearby generator frequencies
figure(3)
plot(TSPAN,OMEGA(:,44),'k',TSPAN,PodOmegar(:,44),'b--',TSPAN,QuadOmegar(44,:),'r--','linewidth',3)
% end of file
