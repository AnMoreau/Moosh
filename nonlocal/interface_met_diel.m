function S=interface(r1,r2,omega)
S=[r1-r2+i*omega,2,2*r2; 2*i*omega*r1,(r2+r1+i*omega),2*i*omega*r2; 2*r1, 2, r2-r1+i*omega]/(r1+r2-i*omega);
end
