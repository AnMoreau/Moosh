

clear all
more off
addpath(genpath(fullfile(fileparts(pwd),'data/')));
addpath('data/');

% Wavelength range

l_min=400;
l_max=800;

% Number of steps
n_steps=101;

lam=linspace(l_min,l_max,n_steps);

% initialization - like Guided modes

for kk=1:n_steps

lambda=lam(kk);
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
% >>>>>>>>>>>>>> First mode retrieval <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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

if (kk==1)
relation=[lambda,modes];
else
if (length(modes)<(size(relation)(2)-1))
modes=[modes,zeros(1,size(relation)(2)-1-length(modes))];
end

if (length(modes)>(size(relation)(2)-1))
relation=[relation,zeros(size(relation)(1),1)];
end

relation=[relation;lambda,modes]
end

end


figure(1)
clf
hold on
for k=1:length(modes)
  color=substr("krgbymcw",rem(k,8)+1,1);
  plot(1.239841974583151e-06./relation(:,1),real(relation(:,k+1)),color,'linewidth',3)
endfor