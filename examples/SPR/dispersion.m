%
% r=dispersion(kx,lambda)
% 
% Computes a function abs(f), whose value is stored in r, so that f(kx,lamba)=0
% is the dispersion relation. 
% Remark: why not having just used the determinant of the scattering matrix ?
% Because it doesn't work very well with the steepest descent, for whatever reason.

function r=dispersion(kx,lambda)

structure

if(length(Type)==3)
Type=[Type(1),Type(2),Type(2:length(Type))];
hauteur=[hauteur(1),hauteur(2)/2,hauteur(2)/2,hauteur(3:length(hauteur))];
end

k0=2*pi/lambda;

if pol==0
  ff=Mu;
else
  ff=Epsilon;
end

g=length(Type);

gamma=zeros(1,g);


% The determination of the square root has to be changed, in order to be able
% to retrieve modes with complex wavevectors.
gamma=sqrt(Epsilon(Type).*Mu(Type).*k0^2-ones(1,g)*kx^2);
gamma(2:g-1)=gamma(2:g-1).*(1-2*(imag(gamma(2:g-1))<0));
gamma(1)=extsqrt(Epsilon(Type(1))*k0.^2-kx.^2);
%gamma(2)=intsqrt(Epsilon(Type(2)).*Mu(Type(2)).*k0^2-kx^2);
gamma(2)=extsqrt(Epsilon(Type(2)).*Mu(Type(2)).*k0^2-kx^2);
%gamma(g-1)=intsqrt(Epsilon(Type(g-1)).*Mu(Type(g-1)).*k0^2-kx^2);
gamma(g-1)=extsqrt(Epsilon(Type(g-1)).*Mu(Type(g-1)).*k0^2-kx^2);
gamma(g)=extsqrt(Epsilon(Type(g))*k0.^2-kx.^2);
% If you experience problems with this program, contact us. A another determination of the square root may have to be considered.

b1=gamma(1)/ff(Type(1));
b2=gamma(2)/ff(Type(2));
T1=[b1+b2,b2-b1;b2-b1,b1+b2]/(2*b2);
for k=1:g-2
  t=exp(i*gamma(k+1)*hauteur(k+1));    
  T{2*k-1}=[0,t;t,0];  
  b1=gamma(k+1)/ff(Type(k+1));
  b2=gamma(k+2)/ff(Type(k+2));  
  T{2*k}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
end
b1=gamma(g-1)/ff(Type(g-1));
b2=gamma(g)/ff(Type(g));
TN=[b1+b2,b2-b1;b2-b1,b1+b2]/(2*b2);

S{1}=cascade(T{1},T{2});
for j=2:size(T,2)-2   
  S{j}= cascade(S{j-1},T{j+1});	
end

% Semi-analytical dispersion relation
r=(-1*TN(1,1)/TN(1,2)-S{size(S,2)}(2,2))*(T1(1,1)/T1(2,1)-S{size(S,2)}(1,1))-S{size(S,2)}(2,1)*S{size(S,2)}(1,2);
r=abs(r);

end

