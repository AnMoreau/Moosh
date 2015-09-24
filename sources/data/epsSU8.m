# La longueur d'onde s'exprime en nanometres
# Wavelength in nanometers

function permittivity=epsSU8(lambda)
l=lambda*1e-3;
%l=lambda;
n=1.566.+0.00796./l.^2.+0.00014./l.^4;
permittivity=n.^2;
endfunction
