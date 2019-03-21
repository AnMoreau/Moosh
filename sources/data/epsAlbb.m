function epsilon=epsAlbb(lambda)

w=6.62606957e-25*299792458/1.602176565e-19/lambda;

f0=0.526;
Gamma0=0.047;
omega_p=14.98;
f=[0.213,0.060,0.182,0.014];
Gamma=[0.312,0.315,1.587,2.145];
omega=[0.163,1.561,1.827,4.495];
sigma=[0.013,0.042,0.256,1.735];

a=sqrt(w*(w+i*Gamma));
a=a.*sign(real(a));
x=(a-omega)./(sqrt(2)*sigma);
y=(a+omega)./(sqrt(2)*sigma);
% Conversion

epsilon=1-omega_p^2*f0/(w*(w+i*Gamma0))+sum(i*sqrt(pi).*f.*omega_p^2./(2*sqrt(2).*a.*sigma).*(faddeeva(x,64)+faddeeva(y,64)));

end
