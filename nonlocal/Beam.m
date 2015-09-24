
clear all;

% Wavelength, nanometers.
lambda=363.8;
%lambda=350;
% Spatial period of the window - discretize means periodize.
%d=15*lambda;
d=7*lambda;
% Waist of the incident gaussian beam.
% Smaller than lambda = almost ponctual source.
w=0.1*lambda;
% Angle of incidence of the beam.
%theta=80*pi/180;
theta=0*pi/180;
% Position of the incident beam with respect to the domain :
% 0.5 means centered on the middle of the domain.
C=0.5;
%C=0.5-200*tan(theta)/d;
%----------- Structure
structure
% Image characteristics :
% Number of pixels horizontally
%nx=lambda*10;
nx=floor(d);
% Number of points for each layer
%ny=[300,repmat([4,9],1,40),4,200];
ny=floor(hauteur);
% Number of modes retained for the description of the field
% so that the last mode has an amplitude < 1e-3
nmod=floor(0.83660*d/w);

%----------- From now on, no need to modify anything---------------

hauteur=hauteur/d;
l=lambda/d;
w=w/d;
% !! Omega_p est aussi normalisée !!
w_p=w_p*d;  
k0=2*pi/l;

En=zeros(sum(ny),nx+1);
En1=zeros(sum(ny),nx+1);
En2=zeros(sum(ny),nx+1);

% Total number of layers
g=length(type);

% Amplitude of the different modes
X=exp(-w^2*pi^2*[-nmod:nmod].^2).*exp(-2*i*pi*[-nmod:nmod]*C);   % amplitude de l'onde incidente
figure(7);
plot(abs(X));
% Scattering matrix corresponding to no interface.
T{1}=[0,1;1,0];

for nm=1:2*nmod+1

% The top medium is assumed to be LOCAL.
  alpha=sqrt(epsilon(type(1)))*k0*sin(theta)+2*pi*(nm-nmod-1);
  gamma(1)=sqrt(epsilon(type(1))*k0^2-alpha^2);

  for k=1:g-1

    if (beta(type(k))==0)
% Layer scattering matrix for a local, dielectric layer
      Kl(k)=0;
      omega(k)=0;
      t=exp(i*gamma(k)*hauteur(k));    
      T{2*k}=[0,t;t,0];
    else
% Layer scattering matrix for a metallic layer.
      Kl(k)=sqrt(alpha^2+((w_p(type(k))^2)*((1/chi_f(type(k)))+(1/(1+chi_b(type(k)))))/beta(type(k))^2));
      omega(k)=alpha^2*(1/(1+chi_f(type(k))+chi_b(type(k)))-1/(1+chi_b(type(k))))/Kl(k);	
      t=exp(i*gamma(k)*hauteur(k));
      l=exp(-Kl(k)*hauteur(k)); 
      T{2*k}=[0 0 t 0; 0 0 0 l; t 0 0 0;0 l 0 0];
    end
      
    gamma(k+1)=sqrt(epsilon(type(k+1))*k0^2-alpha^2);    
% Changing the cut of the square root to achieve perfect stability
    if (imag(gamma(k+1))<0)
      gamma(k+1)=-gamma(k+1);
    end
    if (beta(type(k+1))~=0)
      Kl(k+1)=sqrt(alpha^2+((w_p(type(k+1))^2)*((1/chi_f(type(k+1)))+(1/(1+chi_b(type(k+1)))))/beta(type(k+1))^2));
      omega(k+1)=alpha^2*(1/(1+chi_f(type(k+1))+chi_b(type(k+1)))-1/(1+chi_b(type(k+1))))/Kl(k+1);
    end
	
    if ((beta(type(k))==0)&&(beta(type(k+1))==0))
      b1=gamma(k)/epsilon(type(k));
      b2=gamma(k+1)/epsilon(type(k+1));
      T{2*k+1}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
    else
      if ((beta(type(k))==0)&&(beta(type(k+1))~=0))
% Interface scattering matrix diel -> metal
	b1=gamma(k)/epsilon(type(k));
	b2=gamma(k+1)/epsilon(type(k+1));
	T{2*k+1}=[b1-b2+i*omega(k+1),2*b2,2;2*b1,b2-b1+i*omega(k+1),2; 2*i*omega(k+1)*b1,2*i*omega(k+1)*b2,b1+b2+i*omega(k+1)]/(b1+b2-i*omega(k+1));
      else
% Interface scattering matrix metal -> diel
	b1=gamma(k)/epsilon(type(k));
	b2=gamma(k+1)/epsilon(type(k+1));
	T{2*k+1}=[b1-b2+i*omega(k),-2,2*b2; -2*i*omega(k)*b1,b2+b1+i*omega(k),-2*i*omega(k)*b2; 2*b1,-2, b2-b1+i*omega(k)]/(b1+b2-i*omega(k));
      end
    end
  end    

% Scattering matrix for the last layer
  t=exp(i*gamma(g)*hauteur(g));    
  T{2*g}=[0,t;t,0];

% Once the scattering matrixes have been prepared, now let us combine them.
  
  H{1}=T{2*g};
  A{1}=T{1};
  
  for j=1:size(T,2)-2
    A{j+1}= cascade(A{j},T{j+1}); 	
    H{j+1}= cascade(T{size(T,2)-j},H{j});
  end
  
% And let us compute the intermediate coefficients from the scattering matrixes

  for j=1:size(T,2)-1
    I{j}= intermediaire(A{j},H{size(T,2)-j});
  end

  h=0;
  
  I{2*g}=zeros(2,2);

% Computation of the field in the different layers for one mode (plane wave)

  t=1;
  E=zeros(sum(ny),1);
  for k=1:g
    for m=1:ny(k)
      h=h+hauteur(k)/ny(k);
      	E(t,1)=I{2*k-1}(1,1)*exp(i*gamma(k)*h)+I{2*k}(2+(beta(type(k))~=0),1)*exp(i*gamma(k)*(hauteur(k)-h));
	t=t+1;
      end
      h=0;
  end
  
% For one mode, the image is invariant horizontally
  E=E*exp(i*alpha*[0:nx]./nx);  

% All the images have to be combined using the proper amplitude for each plane wave.
  En=En+X(nm)*E;

end

V=abs(En);

% You have then to visualize the result. May need some optimization, though.

colormap(jet)
imagesc(V);




