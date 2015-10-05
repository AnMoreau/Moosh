%
% [absorb,r,R,t,T]=absorption(theta,lambda,polarization)
%
% This function computes the percentage of the incoming energy
% that is absorbed in each layer when the structure is illuminated
% by a plane wave characterized by an incidence angle "theta" IN RADIANS
% and a wavelength "lambda" IN LENGTH UNIT. It is possible (but not mandatory)
% to specify the polarization (0=s,TE;1=p,TM) of the incoming light.
% If not specified, the value in structure.m will be used.
%
% It returns a vector "absorb" containing the relative absorption in each layer,
% as well as the reflexion and transmission coefficients, as computed by coefficient(),
% for convenience.

function [absorb,r,R,t,T]=absorption(theta,lambda,varargin)

% Calling "structure.m"
structure

if size(varargin)~=0
pol=varargin{1};
end

% In order to get a phase that corresponds to the expected reflected coefficient,
% we make the height of the upper (lossless) medium vanish. It changes only the
% phase of the reflection coefficient, and nothing else (not the absorption, since
% the medium is assumed lossless).
hauteur(1)=0;

if pol==0	%TE
  f=Mu;
else		%TM
  f=Epsilon;
end
k0=2*pi./lambda;
g=length(Type);

alpha=sqrt(Epsilon(Type(1))*Mu(Type(1)))*k0*sin(theta);
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

%first S matrice
T{1}=[0,1;1,0];
for k=1:g-1
% Layer scattering matrix         
	t=exp(i*gamma(k)*hauteur(k));    
	T{2*k}=[0,t;t,0];
% Interface scattering matrix
	b1=gamma(k)/f(Type(k));
	b2=gamma(k+1)/f(Type(k+1));
	T{2*k+1}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
end
t=exp(i*gamma(g)*hauteur(g));    
T{2*g}=[0,t;t,0];

% Once the scattering matrixes have been prepared, now let us combine them
  
  H{1}=T{2*g};
  A{1}=T{1};  
  for j=1:size(T,2)-2
    A{j+1}= cascade(A{j},T{j+1}); 	
    H{j+1}= cascade(T{size(T,2)-j},H{j});
  end  
%reflexion
r=A{size(A,2)}(1,1); 
%transmission
t=A{size(A,2)}(2,1);

% And let us compute the intermediate coefficients from the scattering matrixes
for j=1:size(T,2)-1
  I{j}=[A{j}(2,1),A{j}(2,2)*H{size(T,2)-j}(1,2);A{j}(2,1)*H{size(T,2)-j}(1,1),H{size(T,2)-j}(1,2)]./(1-A{j}(2,2)*H{size(T,2)-j}(1,1));
end
I{2*g}=[I{2*g-1}(1,1)*exp(i*gamma(g)*hauteur(g)),I{2*g-1}(1,2)*exp(i*gamma(g)*hauteur(g));0,0];

%________________________________
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Computation of energy fluxes
%////////////////////////////////
%
% At the top and bottom of each layer
%

w=1;
if (pol == 0)	%TE
  for j=1:2*g		
    poynting(j)=real( (I{j}(1,1)+I{j}(2,1)) * conj( ( I{j}(1,1)-I{j}(2,1) )*gamma(w)/Mu(Type(w)) )*Mu(Type(1))/(gamma(1)));
    w=w+1-mod(j,2);
  end
else		%TM
  for j=1:2*g		
    poynting(j)=real( (I{j}(1,1)-I{j}(2,1)) * conj(I{j}(1,1)+I{j}(2,1))*gamma(w)/Epsilon(Type(w)) *Epsilon(Type(1))/gamma(1));		
    w=w+1-mod(j,2);
  end
end

% Absorption in each layer
tmp=abs(-diff(poynting));
absorb=tmp(1:2:2*g-1);
% reflection coefficient of the whole structure
r=A{size(A,2)}(1,1); 
% transmission coefficient of the whole structure
t=A{size(A,2)}(2,1); 
% Energy reflexion coefficient;
R=abs(r)^2;
% Energy transmission coefficient;
T=abs(t)^2*gamma(g)*f(Type(1))/(gamma(1)*f(Type(g)));

end
