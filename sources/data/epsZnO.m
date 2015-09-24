function eps=epsZnO(lambda)
				%Permittivité de ZnO dans le domaine 375-1100nm
				%Thin Solid Films 510 (2006) 32-38
if ((lambda>350)&&(lambda<1100))
eps=1+2.612*lambda^2/(lambda^2-224.627^2);
else
disp('Out of range for ZnO')
endif

endfunction