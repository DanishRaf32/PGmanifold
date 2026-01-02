function a = symbolic_array(name,n)
% Aufruf: a = symbolic_array(name,n)
% erstellt array der Größe n mit symbolischen Variablen
% falls n eine ganze Zahl ist wird ein Vektor mit symbolischen Variablen
% erstellt

if sum(size(n)) == 2
    a = zeros(n,1);
    syms a;
    for k = 1:n
        string = [name,num2str(k)];
        a(k,1) = sym(string);
    end
else
    a = zeros(n);
    syms a;
    for k = 1:n(1)
        for l = 1:n(2)
            string = [name,num2str(k),num2str(l)];
            a(k,l) = sym(string);
        end
    end
end