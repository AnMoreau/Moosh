function [r,t]=coefficient(theta,lambda)

k0=2*pi/lambda;

structure

g=length(Type);

% Scattering matrix corresponding to no interface.
T{1}=[0,1;1,0];

alpha=sqrt(epsilon(Type(1)))*k0*sin(theta);
gamma(1)=sqrt(epsilon(Type(1))*k0^2-alpha^2);  


for k=1:g-1
    if (Beta(Type(k))==0)
% Layer scattering matrix in the dielectric layer
      Kl(k)=0;
      omega(k)=0;
      t=exp(i*gamma(k)*hauteur(k));    
      T{2*k}=[0,t;t,0];
    else
      Kl(k)=sqrt(alpha^2+((w_p(Type(k))^2)*((1/chi_f(Type(k)))+(1/(1+chi_b(Type(k)))))/Beta(Type(k))^2));
      omega(k)=alpha^2*(1/(1+chi_f(Type(k))+chi_b(Type(k)))-1/(1+chi_b(Type(k))))/Kl(k);	
      t=exp(i*gamma(k)*hauteur(k));
      l=exp(-Kl(k)*hauteur(k)); 
      T{2*k}=[0 0 t 0; 0 0 0 l; t 0 0 0;0 l 0 0];
    end

    gamma(k+1)=sqrt(epsilon(Type(k+1))*k0^2-alpha^2);    

% Changing the cut of the square root to achieve perfect stability
    if (imag(gamma(k+1))<0)
      gamma(k+1)=-gamma(k+1);
    end

% Nonlocal quantities
    if (Beta(Type(k+1))~=0)
      Kl(k+1)=sqrt(alpha^2+((w_p(Type(k+1))^2)*((1/chi_f(Type(k+1)))+(1/(1+chi_b(Type(k+1)))))/Beta(Type(k+1))^2));
      omega(k+1)=alpha^2*(1/(1+chi_f(Type(k+1))+chi_b(Type(k+1)))-1/(1+chi_b(Type(k+1))))/Kl(k+1);
    end
	
    if ((Beta(Type(k))==0)&&(Beta(Type(k+1))==0))
      b1=gamma(k)/epsilon(Type(k));
      b2=gamma(k+1)/epsilon(Type(k+1));
      T{2*k+1}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
    else
      if ((Beta(Type(k))==0)&&(Beta(Type(k+1))~=0))
	% Interface scattering matrix diel -> metal
	b1=gamma(k)/epsilon(Type(k));
	b2=gamma(k+1)/epsilon(Type(k+1));
	T{2*k+1}=[b1-b2+i*omega(k+1),2*b2,2;2*b1,b2-b1+i*omega(k+1),2; 2*i*omega(k+1)*b1,2*i*omega(k+1)*b2, b1+b2+i*omega(k+1)]/(b1+b2-i*omega(k+1));
      else
	% Interface scattering matrix metal -> diel
	b1=gamma(k)/epsilon(Type(k));
	b2=gamma(k+1)/epsilon(Type(k+1));
	T{2*k+1}=[b1-b2+i*omega(k),-2,2*b2;-2*i*omega(k)*b1,b2+b1+i*omega(k),-2*i*omega(k)*b2; 2*b1,-2, b2-b1+i*omega(k)]/(b1+b2-i*omega(k));
      end
    end
end

% Scattering matrix for the last layer
  t=exp(i*gamma(g)*hauteur(g));    
  T{2*g}=[0,t;t,0];

% Once the scattering matrixes have been prepared, now let us combine them.
  A{1}=T{1};
  
for j=1:size(T,2)-2
    A{j+1}= cascade(A{j},T{j+1}); 	
end

 r=A{j+1}(1,1);
% Beware: t^2 is not the energy transmission coefficient - well except for air.
 t=A{j+1}(2,1); 
end