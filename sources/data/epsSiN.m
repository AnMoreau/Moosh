# La longueur d'onde s'exprime en nanometres
# Wavelength in nanometers

function permittivity=epsSiN(lambda)

l=lambda*1e-3;
permittivity=1+2.8939*l^2/(l^2-(139.67e-3)^2);

endfunction
