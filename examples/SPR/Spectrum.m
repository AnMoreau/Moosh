% This program can compute the reflection and transmission coefficients
% as a function of the wavelength, as well as the absorption in a
% specified layer.

clear all
clf
addpath('data/:');

% >>>>>>>>>>>>>>>>>>> Parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Number of the layer where the absorption has to be computed
Layer=3;
% Working incidence angle in degrees
theta=44.6;		
% Polarization - 0 means s or TE; 1 means p or TM.
polarization=1;
% Spectral range in LENGTH UNIT
min=400;
max=800;
% Number of points
Npoints=100;
%_____________________________________________________________________


lambda=700;
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

figure(11)


subplot(2,1,1)
plot(lambda,Ab(:,Layer),'linewidth',2),ylabel('Absorption'),xlabel('Wavelength'),title('Absorption in gold (blue) and chromium (green)');
hold on
plot(lambda,Ab(:,2),'linewidth',2,'g')
hold off


subplot(2,1,2)
plot(lambda,R,'linewidth',2),ylabel('Reflexion'),xlabel('Wavelength'),title('Reflection');


