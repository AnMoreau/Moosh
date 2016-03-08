%_____________________________
% Specific determination of the square root when it has to be changed for the outer media. 
%_____________________________
function truc=intsqrt(z)

truc=sqrt(z);
if (angle(z)<=-pi/5)
  truc=-truc;
end

end
