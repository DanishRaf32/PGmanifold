function J = NumJacobian(f,x_eq,Opts)
%% Parameters
% J: Numerical Jacobi Matrix
% f: function to linearize
% x_eq: Linearization Point

% difference: Difference used to compute the numerical gradient

%% Parse options 
Def.FiniteDifferenceType = 'forward';
Def.FiniteDifferenceStepSize = sqrt(eps);
if exist('Opts','var') && isfield(Opts, 'FiniteDifferenceType') && ...
        strcmp(Opts.FiniteDifferenceType,'central')
    Def.FiniteDifferenceStepSize = eps^(1/3);
end

Def.TypicalX = ones(length(x_eq),1);

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

difference = Opts.FiniteDifferenceStepSize;
type = Opts.FiniteDifferenceType;
TypicalX = Opts.TypicalX;

%% Initialization
% Equilibrium Point
f_eq = f(x_eq);

% Initialize Variables
Jac = zeros(length(f_eq),length(x_eq));
x_minus = x_eq;
x_plus = x_eq;
v = difference*ones(length(x_eq),1);

%% Computation
%Compute the k-th column of the Jacobian
if strcmp(type,'forward')
    delta = v.*signprime(x_eq).*max(abs(x_eq),TypicalX);
    % Forward
    for k = 1:length(x_eq)
        % Generate Numerical Data  
        x_plus(k) = x_plus(k) + delta(k); 
        f_plus = f(x_plus);
    
        %Compute Gradient
        grad = gradient([f_eq f_plus],delta(k));
    
        %Copy Data into the Jacobi Matrix
        Jac(:,k) = grad(:,2);
    
        %Reset Numerical Data
        x_plus(k) = x_plus(k) - delta(k);
    end
elseif strcmp(type,'backward')
    delta = v.*signprime(x_eq).*max(abs(x_eq),TypicalX);
    % Backward
    for k = 1:length(x_eq)
        % Generate Numerical Data
        x_minus(k) = x_minus(k) - delta(k);  
        f_minus = f(x_minus);
    
        %Compute Gradient
        grad = gradient([f_minus f_eq],delta(k));
    
        %Copy Data into the Jacobi Matrix
        Jac(:,k) = grad(:,1);
    
        %Reset Numerical Data
        x_minus(k) = x_minus(k) + delta(k);
    end
elseif strcmp(type,'central')
    delta = v.*max(abs(x_eq),TypicalX);
    % Central
    for k = 1:length(x_eq)
        % Generate Numerical Data
        x_minus(k) = x_minus(k) - delta(k);  
        x_plus(k) = x_plus(k) + delta(k); 
        f_minus = f(x_minus);
        f_plus = f(x_plus);
    
        %Compute Gradient
        grad = gradient([f_minus f_eq f_plus],delta(k));
    
        %Copy Data into the Jacobi Matrix
        Jac(:,k) = grad(:,1);
    
        %Reset Numerical Data
        x_minus(k) = x_minus(k) + delta(k);
        x_plus(k) = x_plus(k) - delta(k);
    end
else
    disp('error: unknown finite difference type')
end

% Return
J = Jac;

    function signp = signprime(x)
        signp = zeros(length(x),1);
        for i = 1:length(x)
            if x(i) == 0
                signp(i) = 1;
            else
                signp(i) = sign(x(i));
            end
        end
    end
end
