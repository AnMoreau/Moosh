% The dispersion relation of the structure can be written
% f(kx,lambda)=0. In order to know

clear all
more off 

addpath("data/:");

% Working wavelength
lambda=700;

k0=2*pi/lambda;
structure

% Discretization of the complex plane
% Number of points along the real axis
nx=300;
% Number of points along the imaginary axis
ny=100;

% Propagation constant window in the complex plane
% Along the real axis
% Effective index range where to look for guided modes
% To find a guided mode if the structure is made only of dielectrics
% use the lowest material index as a lower bound, and the highest index
% as the higher bound.
% With metals, you can choose 0 as a lower bound and there is not perfectly
% correct upper bound when high-k guided modes are studied (gap-plasmons, short-range
% surface plasmons).
xmin=0.5*k0;
xmax=3.2*k0;
% Along the imaginary axis - well it depends on the losses.
ymin=-0.1*k0;
ymax=0.1*k0;

% 
X=linspace(xmin,xmax,nx);
Y=linspace(ymax,ymin,ny);
[ax,ay]=meshgrid(X,Y);
M=ax+i*ay;

for j=1:nx
  j
  for k=1:ny
    r=dispersion(M(k,j),lambda);
% In order to better vizualise where the zeros of the
% above function are, the inverse of its modulus is represented.
    T(k,j)=1/abs(r);
  end
end

figure(1)
hold on
echellex=abs(xmax-xmin)/(nx*k0)   %en k0 par pixel
echelley=abs(ymax-ymin)/(ny*k0)
zero_y=abs(ymin)/(k0*echelley) %en pixel
shading('interp');
colormap(jet(128));
image(min(T,5)/5*128);
hold off


