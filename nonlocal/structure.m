
% Description of the structure :

% 1. Permittivities of the different media considered. 
epsilon=[1,1,epsTiO2(lambda),1.52^2];

% 2. Nonlocal parameters;

%Initialization
chi_f=zeros(1,length(epsilon));
chi_b=chi_f;
w_p=chi_f;
beta=chi_f;

% Actual parameters (gold or silver);
% Medium 2 : silver. Replace Agbb by Aubb if needed for gold parameters.
[epsilon(2),chi_b(2),chi_f(2),w_p(2),beta(2)]=epsAgbb(lambda);

% Purely local solution - put BETA TO ZERO for the material you want
% to be considered as local and uncomment the following line:
%beta=[0,0,0,0];

% 3. How the different media are stacked.

% That is a incident beam in air, on a 10-pattern Bragg mirror, put on top of a metallic layer and on a substrate.
type=[1,repmat([2,3],1,10),2,4];

% 4. Thickness of the different layers

hauteur=[2,repmat([10,42],1,10),10,200];

% End of the structure description
% polarization is always p (TM) otherwise no nonlocal effect takes place.
% Before using Beam.m make sure all the parameters have been changed there.