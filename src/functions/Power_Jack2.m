function [Jf]=Power_Jack2(x)
global data2
n=length(data2.H);
D=data2.D;
H=data2.H;
%A=data.A;
gamma=data2.gamma;
K=data2.K;
omega=data2.omega_R;
nn=length(x);
k=nn/2;
Jf = zeros(nn,nn); 
%f(1:nnodes)=x(nnodes+1:end); %x1dot =x2;
delta = x(1:k);   %x1
cal=zeros(n,n);
for i =1:n
  for j =1:n
     if i==j
         cal(i,j)=0;
     else
         cal(i,j)=K(i,j)*cos(delta(i)-delta(j)-gamma(i,j));
     end
  end
end
fi=sum(cal,2);
Jf(1:k,1:k)= - diag(D./(2*H)) - diag(omega./(2*H).*fi);

end