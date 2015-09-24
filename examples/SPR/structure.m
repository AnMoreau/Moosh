% Characteristics of the different media that are considered.
% Medium #i has for permittivity Epsilon(i) and Mu(i).

Epsilon=[epsbk7(lambda),epsCr(lambda),epsAubb(lambda),epsH2O(lambda),1];
Mu=[1,1,1,1,1];

% How the different media are stacked, beginning with the upper layer
% that is assumed to be lossless for reflection and transmission coefficients.
% Type indicates the number of the medium each layer is made of.
% Layer j has for permittivity Epsilon(Type(j)) and for permeability Mu(Type(j)).
%
% WARNING: The permittivity of the upper layer is assumed to be REAL, 
% which means lossless. If you want the transmission coefficient to have
% any meaning, the last medium has to be lossless too.

Type=[1,2,3,5];

% Height of each layer - in LENGTH UNIT. Be careful to be coherent in all 
% files, and to remember that for most functions providing the permittivity,
% in data, wavelengths are supposed to be in nanometers.
% The height of the first and last layer can be put to zero, except when
% maps of the fields have to be computed.

hauteur=[4000,4,38,1000];

% Polarization considered
% In some functions, this can be overridden.
pol=1;