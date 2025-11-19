%Nonlinear Model Order Reduction of Power GRid Networks using Quadratic
%Manifolds
% Farhana and Danish
% Department of Electrical Engineering, IOT, University of Kashmir, India

% Simulates IEEE 118 Power System Models via MATPOWER and pg_sync_models using
% Euler/ode15s. Reduction via POD and Quadratic ,manifold.
% Srcipt simulates two fault scenerio: 
% 1) A self-clearing 3-phase fault at Bus 1
% 2) A 3-phase fault at bus 1 followed by line-tripping 1-2
% Note: The generator dampings are to be changed to 10pu. 
%   References:
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
global data data2 data3
tspan1 = 0:0.01:5;    % Pre-fault simulation time
tspan2 = 5:0.01:6.5;  % 1500ms Fault time
tspan3 = 6.5:0.01:10; % Post_fault time
%pre-fault 
mpc=case118;    
mpc.ref_freq=60; 
[data,details]=EN_model(mpc); 
%fault
mpc2=case118_3phase_fault; 
mpc2.ref_freq=60;
[data2,details2]=EN_model(mpc2); 
%post-fault
mpc3=case118_3phase_postfault; 
mpc3.ref_freq=60;
[data3,details3]=EN_model(mpc3); 
%System-Parameters
n_oc=length(data.H);   
u=@(t) 1;
solver='Euler'; % Euler, ode15s

%% Pre-Fault
disp('Simulating pre-fault conditions...')
% FOM
f1= @(x) power_func(x); % nonlinear mapping
x01=zeros(2*n_oc,1);
Jack1=@(x) Power_Jack(x);
[yFOM1,xFOM1,FOM_time1,delta1,omega1]=FOM_solve(f1,Jack1,tspan1,x01,solver,n_oc);
% POD
rdefl1=6;
 [Vpod1]=basisRed(xFOM1',rdefl1);
 [yROM_POD1,~,ROM_POD_time1,delta_POD1,omega_POD1]=ROM_solve(f1,Jack1,tspan1,x01,solver,Vpod1,n_oc);
% Quadratic Manifold
lambda1=0.01; % regularization parameter
shat1 = Vpod1'*(xFOM1');
W1 = get_x_sq(shat1')';
err1 = xFOM1'-Vpod1*shat1;
[q,~] = size(xFOM1'); 
[p,~]= size(W1);
Aplus = [W1'; sqrt(lambda1)*eye(p)];
bplus = [err1'; zeros(p,q)];
Vbar1 = (Aplus\bplus)'; %state-dependent basis
[yROM_Quad1,xROM_Quad1,QuadROM1_time,Deltar_Quad1,Omegar_Quad1]=QuadROM_solve(f1,Jack1,tspan1,x01,solver,Vpod1,n_oc,Vbar1,W1);

%%  Fault-Condition at Branch 12
disp('Simulating fault conditions...')
%FOM
f2= @(x) power_func2(x);
Jack2=@(x) Power_Jack2(x);
x02 = xFOM1(end,:)';
[yFOM2,xFOM2,FOM_time2,delta2,omega2]=FOM_solve(f2,Jack2,tspan2,x02,solver,n_oc);
% POD
rdefl2=rdefl1;
[Vpod2]=basisRed(xFOM2',rdefl2); %15
[yROM_POD2,~,ROM_POD_time2,delta_POD2,omega_POD2]=ROM_solve(f2,Jack2,tspan2,x02,solver,Vpod2,n_oc);
% Quadratic Manifold
lambda2=0.01; %1e-1=118, 1=300
shat2 = Vpod2'*(xFOM2');
W2 = get_x_sq(shat2')';
err2 = xFOM2'-Vpod2*shat2;
[q,~] = size(xFOM2'); 
[p,~]= size(W2);
Aplus = [W2'; sqrt(lambda2)*eye(p)];
bplus = [err2'; zeros(p,q)];
Vbar2 = (Aplus\bplus)'; %state-dependent basis
[yROM_Quad2,xROM_Quad2,QuadROM2_time,Deltar_Quad2,Omegar_Quad2]=QuadROM_solve(f2,Jack2,tspan2,x02,solver,Vpod2,n_oc,Vbar2,W2);
% Plotting Outputs
%plotOut(tspan2,yFOM2,yROM_POD2,yROM_Quad2) %check
%% Post-fault scnerio (line 12 cleared)
disp('Simulating post-fault conditions (line 12 tripped)...')
f3= @(x) power_func3(x);
Jack3=@(x) Power_Jack3(x);
x03 = xFOM2(end,:)';
[yFOM3,xFOM3,FOM_time3,delta3,omega3]=FOM_solve(f3,Jack3,tspan3,x03,solver,n_oc);
% POD
rdefl3=rdefl1;
 [Vpod3]=basisRed(xFOM3',rdefl3); 
 [yROM_POD3,~,ROM_POD_time3,delta_POD3,omega_POD3]=ROM_solve(f3,Jack3,tspan3,x03,solver,Vpod3,n_oc);
% Quadratic Manifold
lambda3=0.01; %1e-1=118, 1=300
shat3 = Vpod3'*(xFOM3');
W3 = get_x_sq(shat3')';
err3 = xFOM3'-Vpod3*shat3;
[q,~] = size(xFOM3'); 
[p,~]= size(W3);
Aplus = [W3'; sqrt(lambda3)*eye(p)];
bplus = [err3'; zeros(p,q)];
Vbar3 = (Aplus\bplus)'; %state-dependent basis
[yROM_Quad3,xROM_Quad3,QuadROM3_time,Deltar_Quad3,Omegar_Quad3]=QuadROM_solve(f3,Jack3,tspan3,x03,solver,Vpod3,n_oc,Vbar3,W3);
% Plotting Outputs
%plotOut(tspan3,yFOM3,yROM_POD3,yROM_Quad3) %check
%% combining outputs
YFOM=[yFOM1; yFOM2; yFOM3];
YROM_POD=[yROM_POD1; yROM_POD2; yROM_POD3];
YROM_QUAD=[yROM_Quad1; yROM_Quad2; yROM_Quad3];
TSPAN=[tspan1 tspan2 tspan3];
FOM_TIME=FOM_time1+FOM_time2+FOM_time3;
ROM_POD_TIME = ROM_POD_time1; ROM_POD_time2; ROM_POD_time3;
ROM_QUAD_TIME = QuadROM1_time+QuadROM2_time+QuadROM3_time;
QUAD_TIME = QuadROM1_time + QuadROM2_time + QuadROM3_time;
DELTA=[delta1;delta2;delta3];
DELTA_POD=[delta_POD1;delta_POD2;delta_POD3];
DELTA_QUAD=[Deltar_Quad1, Deltar_Quad2, Deltar_Quad3];
OMEGA=[omega1;omega2;omega3];
OMEGA_POD=[omega_POD1;omega_POD2;omega_POD3];
OMEGA_QUAD=[Omegar_Quad1,Omegar_Quad2,Omegar_Quad3];
% Errors
gen1_POD=abs(OMEGA(:,1)-OMEGA_POD(:,1));
gen5_POD=abs(OMEGA(:,5)-OMEGA_POD(:,5));
gen8_POD=abs(OMEGA(:,8)-OMEGA_POD(:,8));
gen1_QUAD=abs(OMEGA(:,1)-OMEGA_QUAD(1,:)');
gen5_QUAD=abs(OMEGA(:,5)-OMEGA_QUAD(5,:)');
gen8_QUAD=abs(OMEGA(:,8)-OMEGA_QUAD(8,:)');
% norms
norm_gen1POD=norm(gen1_POD)./norm(OMEGA(:,1));
norm_gen5POD=norm(gen5_POD)./norm(OMEGA(:,5))
norm_gen8POD=norm(gen8_POD)./norm(OMEGA(:,8))
norm_gen1QUAD=norm(gen1_QUAD)./norm(OMEGA(:,1))
norm_gen5QUAD=norm(gen5_QUAD)./norm(OMEGA(:,5))
norm_gen8QUAD=norm(gen8_QUAD)./norm(OMEGA(:,8))
figure()
semilogy(TSPAN,gen1_POD,'b',TSPAN,gen5_POD,'r',TSPAN,gen8_POD,'g')
hold on
semilogy(TSPAN,gen1_QUAD,'b--',TSPAN,gen5_QUAD,'r--',TSPAN,gen8_QUAD,'g--')
grid on
xlim([4,10])
xlabel('$t$ [s]','Interpreter','Latex'); ylabel('rel$e_{1}(t)$, rel$e_{5}(t)$, rel$e_{8}(t)$','Interpreter','Latex')
%legend('Linear Manifold','Quadratic Manifold')

%% Plotting Nearby generator frequencies
figure
plot(TSPAN,OMEGA(:,1),'ko',TSPAN,OMEGA_POD(:,1),'b.-',TSPAN,OMEGA_QUAD(1,:),'r--','linewidth',2)
hold on
plot(TSPAN,OMEGA(:,8),'ko',TSPAN,OMEGA_POD(:,8),'b.-',TSPAN,OMEGA_QUAD(8,:),'r--','linewidth',2)
plot(TSPAN,OMEGA(:,5),'ko',TSPAN,OMEGA_POD(:,5),'b.-',TSPAN,OMEGA_QUAD(5,:),'r--','linewidth',2)
xlabel('t(sec)','Interpreter','LaTeX')
ylabel('rad/s','Interpreter','LaTeX')
set(gca,'FontSize',20,'TickLabelInterpreter','Latex')
legend('Reference', 'Linear manifold', 'Quadratic manifold')
% end of file