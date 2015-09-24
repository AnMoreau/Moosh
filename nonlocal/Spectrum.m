% This program can compute the reflection and transmission coefficients
% as a function of the wavelength, as well as the absorption in a
% specified layer.

clear all
%clf
%addpath('data/:');

% >>>>>>>>>>>>>>>>>>> Parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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

lambda=linspace(min,max,Npoints);

for k=1:Npoints		

	[r(k),t(k)]=coefficient(theta*pi/180,lambda(k));

end

% >>>>>>>>>>>>>>>>>>>   Vizualization   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(1)

subplot(3,1,1)
plot(lambda,abs(t),'linewidth',2),ylabel('Absorption'),xlabel('Wavelength'),title('Modulus of the transmission coefficient');

subplot(3,1,2)
plot(lambda,abs(r).^2,'linewidth',2),ylabel('Reflexion'),xlabel('Wavelength'),title('Energy reflection coefficient');

subplot(3,1,3)
plot(lambda,angle(r),'linewidth',2),ylabel('Transmission'),xlabel('Wavelength'),title('Phase of the reflection coefficient');


