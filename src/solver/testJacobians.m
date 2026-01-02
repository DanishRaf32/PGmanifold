clear all; close all; clc;

%% Load nonlinear model
% original order of RC-ladder model
N = 50; 
[f,B,C] = RC_ladder(N);

%% Linearize Model (Symbolic) 
% Define equilibirum at which lit should be linearized
x_eq0 = 0.1*ones(N,1);

% Get the name of the symbolic variable
f_info = functions(f);
str = f_info.function;
f_var = str(3);

%Create Symbolic Variable
sym_var = symbolic_array(f_var,N);

% Derive f wrt sym_var
J = SymJacobian(f,sym_var);

% Substitute
delta_fsym = SubJacobian(J,x_eq0,sym_var.');

%Convert to Standard state space Model
A_linsym = delta_fsym;
B_linsym = B;
C_linsym = C;
D_linsym = zeros(size(C_linsym,1),size(B_linsym,2));
sys_linearizedsym = ss(A_linsym,B_linsym,C_linsym,D_linsym);

%% Linearize Model (Numeric) 
% Derive f at x_eq
Opts.FiniteDifferenceType = 'central';
% Opts.FiniteDifferenceStepSize = 1e-8;
delta_fnum = NumJacobian(f,x_eq0,Opts);
% delta_fnum = NumJacobian(f,x_eq0);

figure;
spy(delta_fnum);

%Convert to Standard state space Model
A_linnum = delta_fnum;
B_linnum = B;
C_linnum = C;
D_linnum = zeros(size(C_linnum,1),size(B_linnum,2));
sys_linearizednum = ss(A_linnum,B_linnum,C_linnum,D_linnum);

%% Simulation 
% Simulation parameters
Ts = 0.1; %Step width
tEnd = 10;
tSim = 0:Ts:tEnd;

%Excitation signal
u = @(t) 0*(t<4) + 1*(t>=4);
% u = @(t) exp(-t);

% f1 = 50; f2 = 1000;
% u = @(t) sin(2*pi*f1*t)+sin(2*pi*f2*t);
% u = @(t) sin(2*pi*f1*t).*sin(2*pi*f2*t);

% Initial state vector
x_0 = 0.1*ones(N,1);

%Simulation of the nonlinear System
f_nonlin = @(t,x) f(x) + B*u(t);
[t_nonlin,x_nonlin] = ode23(f_nonlin,tSim,x_0,struct('RelTol',1e-8));
y_nonlin = C*x_nonlin.';

%Simulation of the symbolically linearized System
f_linsymb = @(t,delta_x) f_nonlin(0,x_eq0) + sys_linearizedsym.A*delta_x + sys_linearizedsym.B*u(t);
[t_symb,delta_x_symb] = ode23(f_linsymb,tSim,x_0 - x_eq0,struct('RelTol',1e-8));
y_symb = sys_linearizedsym.C*(delta_x_symb+(x_eq0*ones(1,length(t_symb))).').';

% Simulation of the numerically linearized System
f_linnum = @(t,delta_x) f_nonlin(0,x_eq0) + sys_linearizednum.A*delta_x+sys_linearizednum.B*u(t);
[t_num,delta_x_num] = ode23(f_linnum,tSim,x_0 - x_eq0,struct('RelTol',1e-8));
y_num = sys_linearizednum.C*(delta_x_num+(x_eq0*ones(1,length(t_num))).').';

%% Plotting
figure(1)
plot(t_nonlin,y_nonlin)
hold on
plot(t_nonlin,y_symb,'k')
plot(t_nonlin,y_num,'r')
legend('Nonlinear System','Symbolic Gradient','Numeric Gradient')
xlabel('t')
ylabel('Output y')
