function M=intermediaire(T,U)

n=min(size(T,1),size(U,1))-1;

m=size(T,1)-n;
p=size(U,1)-n;

A=T(1:m,1:m);
B=T(1:m,m+1:m+n);
C=T(m+1:m+n,1:m);
D=T(m+1:m+n,m+1:m+n);

E=U(1:n,1:n);
F=U(1:n,n+1:n+p);
G=U(n+1:n+p,1:n);
H=U(n+1:n+p,n+1:n+p);

J=inv(eye(n,n)-E*D);
K=inv(eye(n,n)-D*E);


M=[[K*C,K*D*F];[J*E*C,J*F]];

end
