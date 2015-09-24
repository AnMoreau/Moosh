%_____________________________
%définition de la coupure pour les deux milieux externes situés aux extrémités de la structure
%_____________________________
function truc=extsqrt(z)

truc=sqrt(z);
if (angle(z)<=-pi/5)
  truc=-truc;
end

end
