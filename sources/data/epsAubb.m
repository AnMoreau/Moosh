function epsilon=epsAubb(lambda)

w=6.62606957e-25*299792458/1.602176565e-19/lambda;

f0=0.770;
Gamma0=0.050;
omega_p=9.03;
f=[0.054,0.050,0.312,0.719,1.648];
Gamma=[0.074,0.035,0.083,0.125,0.179];
omega=[0.218,2.885,4.069,6.137,27.97];
sigma=[0.742,0.349,0.830,1.246,1.795];

a=sqrt(w*(w+i*Gamma));
a=a.*sign(real(a));
x=(a-omega)./(sqrt(2)*sigma);
#x=x.*sign(real(x));
y=(a+omega)./(sqrt(2)*sigma);
#y=y.*sign(real(x));
% Conversion

#epsilon=1-omega_p^2*f0/(w*(w-i*Gamma0))+sum(i*sqrt(pi).*f*omega_p^2./(2*sqrt(2).*a.*sigma).*(2-cerf(x)-cerf(y)));
epsilon=1-omega_p^2*f0/(w*(w+i*Gamma0))+sum(i*sqrt(pi).*f.*omega_p^2./(2*sqrt(2).*a.*sigma).*(faddeeva(x,64)+faddeeva(y,64)));

end
