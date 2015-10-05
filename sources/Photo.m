% This program computes the short-circuit current for a photovoltaics device.

clear all
addpath(genpath(fullfile(fileparts(pwd),'data/')));
addpath('data/');

% Spectral range -- in nanometers.
% lambda_min>280 nm
lambda_min=375;
% lambda_max < 4000 nm
lambda_max=750;
% Number of points used
% 1000 should be ok for the am1.5 spectrum
n=1000;
% Number of the active layer
active_layer=2;

% Angle of incidence
% Be careful when changing the incidence angle. The incoming power
% is not modified, so...
theta=0;

lambda=linspace(lambda_min,lambda_max,n).';
scc=zeros(n,1);
Ab=scc;
spectrum=scc;

for k=1:n
    absorb=absorption(theta*pi/180,lambda(k));
    scc(k)=solar(lambda(k));
    Ab(k)=absorb(active_layer);
    spectrum(k)=am1_5(lambda(k));
end

max_scc=trapz(lambda,scc)*1e3;
j_sc=trapz(lambda,scc.*Ab)*1e3;
CE=j_sc/max_scc;

subplot(2,1,1)
plot(lambda,spectrum,'linewidth',2)

subplot(2,1,2)
plot(lambda,Ab,'linewidth',2)

disp('Maximum short-circuit current attainable for this spectral range (mA/cm^2)')
disp(max_scc)

disp('Short-circuit current (mA/cm^2)')
disp(j_sc)
disp('Conversion efficiency')
disp(CE)

