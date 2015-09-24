% This program tries to find all the guided modes of a given structure. It must be provided
% with a effective index range, though (ie a range of propagation constants where to look
% for solutions of the dipersion relation). 
%
% Once the modes are found, try "Profile(modes(k))" to get the profile of the k-th mode. 

clear all
addpath("data/:");
% Working wavelength
lambda=700;
k0=2*pi/lambda;

% Effective index range where to look for guided modes
% To find a guided mode if the structure is made only of dielectrics
% use the lowest material index as a lower bound, and the highest index
% as the higher bound.
% With metals, you can choose 0 as a lower bound and there is not perfectly
% correct upper bound when high-k guided modes are studied (gap-plasmons, short-range
% surface plasmons).

n_min=1.05;
n_max=2.352;

% >>>>>>>>>>>>>> Advanced parameters <<<<<<<<<<<<<<<<<<<<<<<<
% These parameters may need to be changed if the program can't
% find satisfying solutions. No rule here, you just have to guess
% what went wrong and adapt the parameters as a consequence.
% (maybe the step is too small and the descent can't reach a solution in
% 1000 steps, which is the upper limit / maybe precision is not small
% enough so that the program finds way too many solutions).
% Number of times the steepest descent has to be tried
% Nstart has to be much larger than the number of modes
Nstart=50;
% Initial step for the steepest descent - increase it if the descent
% is too slow or if it doesn't reach a solution in less than 1000 steps
% (in which case it returns 1000 in "msg", meaning no solution was found).
step=0.05*k0;
% Value below which the function "dispersion" should be considered null.
% If it is not low enough, the values of the propagation constants may change
% when "precision" is changed.
precision=1e-13;

% >>>>>>>>>>>>>> Mode retrieval <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

kx_start=linspace(n_min*k0,n_max*k0,Nstart);
structure
% Using the semi-analytical dispersion relation, a few layers are required.
if(length(Type)<3)
disp('You have to put at least 3 layers in the structure to find a guided mode');
else
if(length(Type)==3)
Type=[Type(1),Type(2),Type(2:length(Type))];
hauteur=[hauteur(1),hauteur(2)/2,hauteur(2)/2,hauteur(3:length(hauteur))];
end

modes=[0];
for k=1:length(kx_start)
% Steepest descent in the complex plane to find a solution of the dispersion relation
  [kx,steps]=descent(lambda,kx_start(k),step,precision,@dispersion);
%  steps
% A guided mode has been found if the solution is far enough from any solution and 
% if the steepest descent has succeed (ie not returned 1000)
  if ((steps~=1000)&&(min(abs(modes-kx))>2e-4*k0))
   modes=[modes,kx];
   disp('New mode found -- effective index')
   tmp=kx/k0;
   disp(tmp)
  end
end
modes=modes(2:length(modes));

disp('Effective indexes of the modes')
disp(modes.'/k0)
disp('Propagation constants stored in variable "modes"')

end
