function [ df ] = SubJacobian(J,x_eq,sym_var)
%sym_varstring = mat2str(sym_var);
for k = 1:length(sym_var)
    %subs(sym_var(k),x_eq(k))
    sym_var(k) = x_eq(k);
end
%sym_var = x_eq;
J_sub = subs(J,sym_var(1));
for k = 2:length(sym_var)
    J_sub = subs(J_sub,sym_var(k));
end
df = double(J_sub);
end
