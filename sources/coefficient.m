%
% [r,R,t,T]=coefficient(theta,lambda,polarization)
%
% This function computes the reflection and transmission coefficients
% of the structure described in structure.m when it is illuminated with
% a plane wave (incidence angle theta in RADIANS, wavelength in vacuum
% lambda in LENGTH UNIT), and the polarization (eventually, otherwise the
% polarization specified in structure.m is used.
% r and t are the amplitude coefficients (complex quantities)
% R and T are the energy coefficients (real quantities)
%
% The transmission coefficients have a meaning only if the lower medium
% is lossless, or they have not meaning AT ALL.

function [r,R,t,T]=coefficient(theta,lambda,varargin)

structure

if size(varargin)!=0
 pol=varargin{1};
end

% In order to get a phase that corresponds to the expected reflected coefficient,
% we make the height of the upper (lossless) medium vanish. It changes only the
% phase of the reflection coefficient.
hauteur(1)=0;

if pol==0	%TE
  f=Mu;
else		%TM
  f=Epsilon;
end

k0=2*pi./lambda;
g=length(Type);

alpha= sqrt(Epsilon(Type(1))*Mu(Type(1)))*k0*sin(theta);
gamma=sqrt(Epsilon(Type).*Mu(Type).*k0^2-ones(1,g)*alpha^2);

% Be cautious if the upper medium is a negative index one.
if (real(Epsilon(Type(1)))<0) && (real(Mu(Type(1)))<0 )
    gamma(1)=-gamma(1);
end
% Changing the determination of the square root to achieve perfect stability
if (g>2)
  gamma(2:g-1)=gamma(2:g-1).*(1-2*(imag(gamma(2:g-1))<0));
end
% Outgoing wave condition for the last medium
if (real(Epsilon(Type(g)))<0) && (real(Mu(Type(g)))<0) && (real(sqrt(Epsilon(Type(g))*Mu(Type(g))*k0^2-alpha^2))~=0) 
  gamma(g)=-sqrt(Epsilon(Type(g))*Mu(Type(g))*k0^2-alpha^2);
else
  gamma(g)=sqrt(Epsilon(Type(g))*Mu(Type(g))*k0^2-alpha^2);
end

T{1}=[0,1;1,0];
%S matrices computation
for k=1:g-1
% Layer scattering matrix 
  t=exp(i*gamma(k)*hauteur(k));    
  T{2*k}=[0,t;t,0];
% Interface scattering matrix
  b1=gamma(k)/f(Type(k));
  b2=gamma(k+1)/f(Type(k+1));
  T{2*k+1}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);  
end
% Scattering matrix for the last layer  
t=exp(i*gamma(g)*hauteur(g));    
T{2*g}=[0,t;t,0];
% Once the scattering matrixes have been prepared, now let us combine them
  
A{1}=T{1};  
for j=1:size(T,2)-2
  A{j+1}= cascade(A{j},T{j+1}); 	
end
  
% reflection coefficient of the whole structure
r=A{size(A,2)}(1,1); 
% transmission coefficient of the whole structure
t=A{size(A,2)}(2,1); 
% Energy reflexion coefficient;
R=abs(r)^2;
% Energy transmission coefficient;
T=abs(t)^2*gamma(g)*f(Type(1))/(gamma(1)*f(Type(g)));

end
