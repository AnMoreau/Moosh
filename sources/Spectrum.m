% This program can compute the reflection and transmission coefficients
% as a function of the wavelength, as well as the absorption in a
% specified layer.

clear all
clf
addpath(genpath(fullfile(fileparts(pwd),'data/')));
addpath('data/')


% >>>>>>>>>>>>>>>>>>> Parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Number of the layer where the absorption has to be computed
Layer=2;
% Working incidence angle in degrees
theta=45;		
% Polarization - 0 means s or TE; 1 means p or TM.
polarization=0;
% Spectral range in LENGTH UNIT
min=375;
max=750;
% Number of points
Npoints=100;
%_____________________________________________________________________


% Structure geometrical parameters
structure

Ab=zeros(Npoints,length(Type));
lambda=linspace(min,max,Npoints);

for k=1:Npoints		

% Call to absorption
	[absorb,r(k),R(k),t(k),T(k)]=absorption(theta*pi/180,lambda(k),polarization);
	Ab(k,:)=absorb;
% If only the reflection coefficient is important, use coefficient() instead
%	[r(k),R(k),t(k),T(k)]=coefficient(theta*pi/180,lambda(k),polarization);

end

% >>>>>>>>>>>>>>>>>>>   Vizualization   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(1)

subplot(3,1,1)
gg=sprintf('Absorption in layer %d',Layer);
plot(lambda,Ab(:,Layer),'linewidth',2),ylabel('Absorption'),xlabel('Wavelength'),title(gg);

subplot(3,1,2)
plot(lambda,R,'linewidth',2),ylabel('Reflexion'),xlabel('Wavelength'),title('Reflection');

subplot(3,1,3)
plot(lambda,T,'linewidth',2),ylabel('Transmission'),xlabel('Wavelength'),title('Phase of the reflection coefficient');

% For test reasons - R+T+Absorption = 1
% R+T+Ab(:,Layer).'-1

