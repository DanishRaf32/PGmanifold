function [f]=power_func3(x)
global data3
n=length(data3.H);
D=data3.D;
H=data3.H;
A=data3.A;
gamma=data3.gamma;
K=data3.K;
omega=data3.omega_R;
nnodes=numel(x)/2;      
f = zeros(size(x)); 
f(1:nnodes)=x(nnodes+1:end); %x1dot =x2;
delta = x(1:nnodes);   %x1
cal=zeros(n,n);
for i =1:n
  for j =1:n
     if i==j
         cal(i,j)=0;
     else
         cal(i,j)=K(i,j)*sin(delta(i)-delta(j)-gamma(i,j));
     end
  end
end
fi=sum(cal,2);
f(nnodes+1:end)= - D./(2*H).*f(1:nnodes) + omega./(2*H).*A -...
                 omega./(2*H).*fi;

end