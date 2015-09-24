% This program can compute the reflection and transmission coefficients
% as a function of the incidence angle, as well as the absorption in a
% specified layer.

clear all

%addpath('data/:');

% >>>>>>>>>>>>>>>>>>> Parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Number of the layer where the absorption has to be computed
Layer=2;
% Workgin wavelength
lambda=500;		
% Polarization - 0 means s or TE; 1 means p or TM.
polarization=0;
% Angular range in degrees here
min=0;
max=89;
% Number of points
Npoints=200;
%_____________________________________________________________________


% Structure geometrical parameters
structure

Ab=zeros(Npoints,length(type));
theta=linspace(min,max,Npoints);

for k=1:Npoints		

	[r(k),t(k)]=coefficient(theta(k)*pi/180,lambda,polarization);

end

% >>>>>>>>>>>>>>>>>>>   Vizualization   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(1)

subplot(3,1,1)
gg=sprintf('Modulus of the transmission coefficient',Layer);
plot(theta,abs(t),'linewidth',2),ylabel('Absorption'),xlabel('Angle (degrees)'),title(gg);


subplot(3,1,2)
plot(theta,abs(r).^2,'linewidth',2),ylabel('Coefficient'),xlabel('Angle'),title('Energy reflection coefficient');

subplot(3,1,3)
plot(theta,angle(r),'linewidth',2),ylabel('Phase (in radians)'),xlabel('Angle in degree'),title('Phase of the reflection coefficient');

% For test reasons - R+T+Absorption = 1
% R+T+Ab(:,Layer).'-1
