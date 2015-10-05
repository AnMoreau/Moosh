function epsilon=epsAgbb(lambda)

w=6.62606957e-25*299792458/1.602176565e-19/lambda;

f0=0.821;
Gamma0=0.049;
omega_p=9.01;
f=[0.050,0.133,0.051,0.467,4.000];
Gamma=[0.189,0.067,0.019,0.117,0.052];
omega=[2.025,5.185,4.343,9.809,18.56];
sigma=[1.894,0.665,0.189,1.170,0.516];
#f=[0.050,0.133,0.051,0.467];
#Gamma=[0.189,0.067,0.019,0.117];
#omega=[2.025,5.185,4.343,9.809];
#sigma=[1.894,0.665,0.189,1.170];

a=sqrt(w*(w+i*Gamma));
a=a.*sign(real(a));
x=(a-omega)./(sqrt(2)*sigma);
y=(a+omega)./(sqrt(2)*sigma);

% Conversion
aha=i*sqrt(pi).*f.*omega_p^2./(2*sqrt(2).*a.*sigma).*(faddeeva(x,64)+faddeeva(y,64));
epsilon=1-omega_p^2*f0/(w*(w+i*Gamma0))+sum(aha);
#(2-cerf(x)-cerf(y)));

end
