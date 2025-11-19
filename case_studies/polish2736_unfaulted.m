%Nonlinear Model Order Reduction of Power GRid Networks using Quadratic
%Manifolds
% Farhana and Danish
% Department of Electrical Engineering, IOT, University of Kashmir, India

% This script is used for nonlinear MOR of Polish 2736 Bus System via linear 
% and quadratic manifolds

% References

% [1]Nishikawa, T. and Motter, A.E., 2015. Comparative analysis of existing models for power-grid synchronization. New Journal of Physics, 17(1), p.015012.
% [2]Geelen, R., Wright, S. and Willcox, K., 2023. Operator inference for non-intrusive model reduction with quadratic manifolds. Computer Methods in Applied Mechanics and Engineering, 403, p.115717. 
%%
clear;clc;close all
mpc=case2736sp;        %57, 118, 145, 300, 1888, 2736sp, test_system_10gen (see data file)
mpc.ref_freq=60;  % reference frequency 
global data data2 
%% Full Order MOdel
[data,details]=EN_model(mpc); %EN, SM (H,D,A,K,gamma,omega_R)
n_oc=length(data.H); %No of oscillators (FOM Size =2*n_oc)
dt=0.001;tf=5;
tspan=0:dt:tf;
f= @(x) power_func(x);
Jack=@(x) Power_Jack(x); %Jacobian function
x0=zeros(2*n_oc,1); % initial conditions
solver='Euler'; %ode45,Euler, ode15s
disp('Solving FOM...')
[yFOM,xFOM,FOM_time,delta_pre,omega_pre]=FOM_solve(f,Jack,tspan,x0,solver,n_oc);
disp('Done!')
%% Simulation of ROMs for diffrerent dimensions
lambda=0.001; % regularization parameter
rmax=90;  % max basis size
E_lin=zeros(1,rmax); % energy captured via linear basis
E_nonlin=zeros(1,rmax); % energy captured via quad basis
err_lin=zeros(2,rmax); % err in linear basis
err_nonlin=zeros(2,rmax); % err in quad basis
CPU_time=zeros(2,rmax); % err in quad basis
yROMs=cell(1,rmax); % store ROM from linear manifold
xROMs=cell(1,rmax);
yQuadROMs=cell(1,rmax); % store ROM from quad manifold
xQuadROMs=cell(1,rmax);
Vpods=cell(1,rmax); Vbars=cell(1,rmax);
rec_lin=cell(1,rmax);rec_nonlin=cell(1,rmax);
i=1;
[U,S,V]=svd(xFOM'); % SVD
rdim=1:1:rmax; % select ROM dimensions
ww = waitbar(0, 'Simulating ROMs...');  % Initialize waitbar
%main loop
for r=rdim
    waitbar(i/length(rdim), ww, sprintf('Solving ROMs: %d%%', round(100*(i/length(rdim)))));  % Update waitbar
    %disp(['Solving ROMs : ', num2str(i),' of ', num2str(length(rdim))])
    Vpod=U(:,1:r); % linear basis truncation; 
    Vpods{i}=Vpod;
    [yROM,xROM,ROM_time,deltar,omegar]=ROM_solve(f,Jack,tspan,x0,solver,Vpod,n_oc);
    yROMs{i}=yROM;xROMs{i}=xROM;
    shat = Vpod'*(xFOM'); %linear basis
    num_lin=Vpod*shat; %reconstructed linear basis
    rec_lin{i}=num_lin;
    E_lin(1,i)=norm(num_lin,'fro')^2 /norm(xFOM','fro')^2; %energy linear basis
    err_lin(1,i)=norm(abs(yFOM(:,1)-yROM(:,1)))./norm(yFOM(:,1)); %delta
    err_lin(2,i)=norm(abs(yFOM(:,2)-yROM(:,2)))./norm(yFOM(:,2)); %omega
    W = get_x_sq(shat')';
    err = xFOM'-Vpod*shat;
    [q,~] = size(xFOM'); 
    [p,~]= size(W);
    Aplus = [W'; sqrt(lambda)*eye(p)];
    bplus = [err'; zeros(p,q)];
    Vbar = (Aplus\bplus)'; %state-dependent basis
    Vbars{i}=Vbar;
    [yQuadROM,xQuadROM,QuadROM_time,Quaddeltar,Quadomegar]=QuadROM_solve(f,Jack,tspan,x0,solver,Vpod,n_oc,Vbar,W);
    yQuadROMs{i}=yQuadROM;xQuadROMs{i}=xQuadROM;
    num_nonlin=num_lin + Vbar*(W); rec_nonlin{i}=num_nonlin;
    E_nonlin(1,i)=norm(num_nonlin,'fro')^2 /norm(xFOM','fro')^2;
    err_nonlin(1,i)=norm(abs(yFOM(:,1)-yQuadROM(:,1)))./norm(yFOM(:,1)); %delta
    err_nonlin(2,i)=norm(abs(yFOM(:,2)-yQuadROM(:,2)))./norm(yFOM(:,2)); %omega
    CPU_time(1,i)=ROM_time;CPU_time(2,i)=QuadROM_time;
    i=i+1;
end

%% Plotting
pp=30; % select max rank to plotting [1 to rmax]
r=3 % reduced basis

;%6=118, 4=300, 3=2736  % select fixed rank for plotting outputs
plot118(tspan,rdim,E_lin,E_nonlin,err_lin,err_nonlin,CPU_time,yFOM,yROMs,yQuadROMs,xFOM,xROMs,xQuadROMs,pp,r)
err_deltar=norm(abs(yFOM(:,1) - yROMs{1,3}(:,1))); 
err_omegar=norm(abs(yFOM(:,2) - yROMs{1,3}(:,2))); 
err_Quaddeltar=norm(abs(yFOM(:,1) - yQuadROMs{1,3}(:,1))); 
err_Quadomegar=norm(abs(yFOM(:,2) - yQuadROMs{1,3}(:,2))); 
% end of file