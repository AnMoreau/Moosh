%_____________________________
%définition de la coupure pour les deux milieux situés àl'intérieur de la structure entre les deux milieux externes
%_____________________________
function truc=intsqrt(z)

truc=sqrt(z);
if (angle(z)<=-pi/5)
  truc=-truc;
end

end
