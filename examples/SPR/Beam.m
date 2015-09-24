

clear all
figure(13)

addpath("data/:");

% Wavelength in vacuum.

lambda=600;

% Structure and electromagnetic parameters at the working wavelength
structure

% Spatial window size
d=70*lambda;

% incident beam width
w=10*lambda;

% Angle of incidence - in radians.
theta=44.6*pi/180;

% Characteristics of the image
% Position of the beam in the spatial window
% 0.5 means centered.
C=0.4;

% Number of pixels horizontaly
nx=floor(d/3);

% Number of pixels verticaly
ny=floor(hauteur);

% Number of modes retained for the description of the field
% so that the last mode has an amplitude < 1e-3 - you may want
% to change it if the structure present reflexion coefficients
% that are subject to very swift changes with the angle of incidence.

nmod=floor(0.83660*d/w);

%----------- Do not touch this part ---------------
% But have a look at the vizualisation part at the very end

hauteur=hauteur/d;
l=lambda/d;
w=w/d;

if pol==0
  f=Mu;
else
  f=Epsilon;
end
  
k0=2*pi/l;

En=zeros(sum(ny),nx+1);
%En1=zeros(sum(ny),nx+1);
%En2=zeros(sum(ny),nx+1);

% Total number of layers
g=length(Type);

% Amplitude of the different modes
X=exp(-w^2*pi^2*[-nmod:nmod].^2).*exp(-2*i*pi*[-nmod:nmod]*C);   % amplitude de l'onde incidente

% Scattering matrix corresponding to no interface.
T{1}=[0,1;1,0];
for nm=1:2*nmod+1
alpha= sqrt(Epsilon(Type(1))*Mu(Type(1)))*k0*sin(theta)+2*pi*(nm-nmod-1);
gamma=sqrt(Epsilon(Type).*Mu(Type).*k0^2-ones(1,g)*alpha^2);

% Taking the ingoing wave condition into account for the above medium.
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
% Once the scattering matrixes have been prepared, now let us combine them.
  
  H{1}=T{2*g};
  A{1}=T{1};
  
  for j=1:size(T,2)-2
    A{j+1}= cascade(A{j},T{j+1}); 	
    H{j+1}= cascade(T{size(T,2)-j},H{j});
  end
  
% And let us compute the intermediate coefficients from the scattering matrixes

  for j=1:size(T,2)-1
	%I{j}= intermediaire(A{j},H{size(T,2)-j});
	I{j}=[A{j}(2,1),A{j}(2,2)*H{size(T,2)-j}(1,2);A{j}(2,1)*H{size(T,2)-j}(1,1),H{size(T,2)-j}(1,2)]./(1-A{j}(2,2)*H{size(T,2)-j}(1,1));
  end

  h=0;
  
  I{2*g}=zeros(2,2);

% Computation of the field in the different layers for one mode (plane wave)

  t=1;
  E=zeros(sum(ny),1);
  for k=1:g
    for m=1:ny(k)
      h=h+hauteur(k)/ny(k);
      	E(t,1)=I{2*k-1}(1,1)*exp(i*gamma(k)*h)+I{2*k}(2,1)*exp(i*gamma(k)*(hauteur(k)-h));
	t=t+1;
	end
    h=0;
  end
  
% For one mode, the image is invariant horizontally
  E=E*exp(i*alpha*[0:nx]./nx);  

% All the images have to be combined using the proper amplitude of each plane wave.
  En=En+X(nm)*E;

end

% -------------- Vizualisation --------------------------
% Can be modified at will.
%

figure(1)
%V=abs(En);
V=real(En);
colormap(jet);

% Generating an image with octave

imwrite(uint8(V/abs(max(max(V)))*63)+1,jet,"spr.jpg");

imagesc(V/max(max(V))*128);
