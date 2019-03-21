function epsilon=epsNibb(lambda)

w=6.62606957e-25*299792458/1.602176565e-19/lambda;

f0=0.083;
Gamma0=0.022;
omega_p=15.92;
f=[0.357,0.039,0.127,0.654];
Gamma=[2.820,0.120,1.822,6.637];
omega=[0.317,1.059,4.583,8.825];
sigma=[0.606,1.454,0.379,0.510];

a=sqrt(w*(w+i*Gamma));
a=a.*sign(real(a));
x=(a-omega)./(sqrt(2)*sigma);
y=(a+omega)./(sqrt(2)*sigma);
% Conversion

epsilon=1-omega_p^2*f0/(w*(w+i*Gamma0))+sum(i*sqrt(pi).*f.*omega_p^2./(2*sqrt(2).*a.*sigma).*(faddeeva(x,64)+faddeeva(y,64)));

end
