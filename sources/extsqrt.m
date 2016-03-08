%_____________________________
% Specific determination of the square root for inner media. 
%_____________________________
function truc=extsqrt(z)

truc=sqrt(z);
if (angle(z)<=-pi/5)
  truc=-truc;
end

end
