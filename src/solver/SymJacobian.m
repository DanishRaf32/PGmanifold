function [ J ] = SymJacobian(f,sym_var)
g = f(sym_var);
J =  jacobian(g,sym_var);
end

