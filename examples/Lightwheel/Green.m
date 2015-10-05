% Computes the field emitted by an oscillating current in a infinitely thin wire
% perpendicular to the plane.

clear all

% Wavelength in nm
lambda=300;

structure

%number of the layer where the source is located
source=3;		

% Position of the source inside the layer
Ys=hauteur(source)/2;

%Horizontal position of the source in the selected layer 
% (0.5 = midle)

Xs=0.5;		

% Spatial window size
% Usually a good idea if it is actually not a multiple of lambda
d=40.56*lambda;

% Number of modes retained 
% Each mode can be related to a given spatial period. Be sure that
% there are enough modes to reach the spatial periods corresponding
% to physically meaningful modes of the structure. 

nmod=100;

%----------- Do not touch this part ---------------
Type=[Type(1:source),Type(source:length(Type))];
hauteur=[hauteur(1:source),hauteur(source:length(hauteur))];
hauteur(source)=Ys;
hauteur(source+1)=hauteur(source+1)-Ys;
%---------------------------------------------------------

% Number of pixels horizontaly
nx=d/8;

% Number of pixels verticaly
ny=floor(hauteur./4);


%----------- Do not touch this part ---------------

hauteur=hauteur/d;
l=lambda/d;

if pol==0
  f=Mu;
else
  f=Epsilon;
end
  
k0=2*pi/l;
N=length(Type);

En=zeros(sum(ny),nx+1);

for nm=1:2*nmod+1

alpha=2*pi*(nm-nmod-1);
gamma=sqrt(Epsilon(Type).*Mu(Type).*k0^2-ones(1,N)*alpha^2);

if (N>2)
  gamma(2:N-1)=gamma(2:N-1).*(1-2*(imag(gamma(2:N-1))<0));
end

% Outgoing wave condition for the last medium
    if (real(Epsilon(Type(N)))<0) && (real(Mu(Type(N)))<0) && (real(sqrt(Epsilon(Type(N))*Mu(Type(N))*k0^2-alpha^2))~=0) 
	gamma(N)=-sqrt(Epsilon(Type(N))*Mu(Type(N))*k0^2-alpha^2);
    else
      gamma(N)=sqrt(Epsilon(Type(N))*Mu(Type(N))*k0^2-alpha^2);
    end


n1=source-1;
n2=N-n1-2;

%Calculation above the source 

  for k=1:n1
    t=exp(i*gamma(k)*hauteur(k));    
    T1{2*k-1}=[0,t;t,0];
    b1=gamma(k)/f(Type(k));
    b2=gamma(k+1)/f(Type(k+1));
    T1{2*k}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
    end
    t=exp(i*gamma(n1+1)*hauteur(n1+1));    
    T1{2*n1+1}=[0,t;t,0];

% Once the scattering matrixes have been prepared, now let us combine them.
  
  H1{1}=T1{size(T1,2)};
  A1{1}=T1{1};
  
  for j=1:size(T1,2)-1
    A1{j+1}= cascade(A1{j},T1{j+1}); 	
    H1{j+1}= cascade(T1{size(T1,2)-j},H1{j});
  end
  
% And let us compute the intermediate coefficients from the scattering matrixes

  clear I1
  for j=1:size(T1,2)-1
	  I1{j}=[A1{j}(2,1)*H1{size(T1,2)-j}(1,1),H1{size(T1,2)-j}(1,2);A1{j}(2,1),A1{j}(2,2)*H1{size(T1,2)-j}(1,2)]./(1-A1{j}(2,2)*H1{size(T1,2)-j}(1,1));
  end

  r1=A1{length(A1)}(2,2);

% Computation of the field in the different layers for one mode (plane wave)

I1=[[0,0;0,0],I1,[0,1;0,r1*exp(2*i*gamma(n1+1)*hauteur(n1+1))]];

% Scattering matrix for the last layer
    t=exp(i*gamma(n1+2)*hauteur(n1+2));    
    T2{1}=[0,t;t,0];

%______________________________________________________________________________    
for k=1:n2 %attention ...
    p=n1+k+1;
    b1=gamma(p)/f(Type(p));
    b2=gamma(p+1)/f(Type(p+1));
    T2{2*k}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
    t=exp(i*gamma(p+1)*hauteur(p+1));    
    T2{2*k+1}=[0,t;t,0];		
end
  
  H2{1}=T2{length(T2)};
  A2{1}=T2{1};
  
  for j=1:size(T2,2)-2
    A2{j+1}= cascade(A2{j},T2{j+1}); 	
    H2{j+1}= cascade(T2{size(T2,2)-j},H2{j});
  end
  
% And let us compute the intermediate coefficients from the scattering matrices
  clear I2
  for j=1:size(T2,2)-1
    I2{j}=[A2{j}(2,1)*H2{size(T2,2)-j}(1,1),H2{size(T2,2)-j}(1,2);A2{j}(2,1),A2{j}(2,2)*H2{size(T2,2)-j}(1,2)]./(1-A2{j}(2,2)*H2{size(T2,2)-j}(1,1));
  end
  I2=[[r1*exp(i*2*hauteur(n1+1)*gamma(n1+1)),0;1,0],I2,[0,0;0,0]];

  h=0;
  
r2=A2{length(A2)}(1,1);


% In the layer of the source

ex=i*exp(i*alpha*Xs);
U1=-ex/(gamma(n1+1)*(r1-1)+gamma(n1+2)*(r2-1)*(1+r1)/(1+r2));
D2=-ex/(gamma(n1+1)*(r1-1)*(1+r2)/(1+r1)+gamma(n1+2)*(r2-1));


% Calculation of the field above the source <<<<<<<<<<<<<<<<<
%if (nm>1)
%clear E
%clear E1
%clear E2
%end

  t1=1;
  E=zeros(sum(ny),1);
E1=zeros(sum(ny(1:n1+1)),1);
E2=zeros(sum(ny(n1+2:N)),1);
  for k=1:n1+1
    for m=1:ny(k)
      h=h+hauteur(k)/ny(k);
      	E1(t1,1)=I1{2*k}(1,2)*U1*exp(i*gamma(k)*(hauteur(k)-h))+I1{2*k-1}(2,2)*U1*exp(i*gamma(k)*(h));
	t1=t1+1;
	end
    h=0;
  end

%calculation of the field below the source <<<<<<<<<<<<<<<<<
  t2=1;
  for k=1:n2+1
  p=n1+1+k;
    for m=1:ny(p)
      h=h+hauteur(p)/ny(p);
      	E2(t2,1)=I2{2*k}(1,1)*D2*exp(i*gamma(p)*(hauteur(p)-h))+I2{2*k-1}(2,1)*D2*exp(i*gamma(p)*h);
	t2=t2+1;
	end
    h=0;
  end


 % All the images have to be combined using the proper amplitude of each plane wave.

tmp=exp(i*alpha*[0:nx]./nx);
En=En+[E1*tmp;E2*tmp];

end

% -------------- Visualization --------------------------

figure(1)
V=abs(En).^2;
colormap(jet);

% Generating an image with octave
% imwrite(uint8(V/abs(max(max(V)))*63)+1,jet,'light_wheel.jpg');

%image_viewer("eog %s")
image(V/max(max(V))*128);
colorbar;
