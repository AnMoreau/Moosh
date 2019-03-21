function epsilon=epsPtbb(lambda)

w=6.62606957e-25*299792458/1.602176565e-19/lambda;

f0=0.333;
Gamma0=0.080;
omega_p=9.59;
f=[0.186,0.665,0.551,2.214];
Gamma=[0.498,1.851,2.604,2.891];
omega=[0.782,1.317,3.189,8.236];
sigma=[0.031,0.096,0.766,1.146];

a=sqrt(w*(w+i*Gamma));
a=a.*sign(real(a));
x=(a-omega)./(sqrt(2)*sigma);
y=(a+omega)./(sqrt(2)*sigma);
% Conversion

epsilon=1-omega_p^2*f0/(w*(w+i*Gamma0))+sum(i*sqrt(pi).*f.*omega_p^2./(2*sqrt(2).*a.*sigma).*(faddeeva(x,64)+faddeeva(y,64)));

end
