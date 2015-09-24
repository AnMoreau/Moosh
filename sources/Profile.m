% 
%     Profile(kx,lambda)
% 
% This function draws the profile of a particular guided mode, that has a wavevector kx (complex)
% for a given wavelength in vaccuum (real, in LENGTH UNIT). Guidedmodes.m has to run first, then
% you can call Profile that way: Profile(modes(n),lambda) where lambda is replaced by the 
% wavelength's value and n is the number of the mode you want to study.

function Profile(kx,lambda)

addpath("data/:");

structure

if(length(Type)<3)
disp('You have to put at least 3 layers in structure to find a guided mode');
else

if(length(Type)==3)
Type=[Type(1),Type(2),Type(2:length(Type))];
hauteur=[hauteur(1),hauteur(2)/2,hauteur(2)/2,hauteur(3:length(hauteur))];
end

% Number of points for the representation of the field.
ny=floor(hauteur.*4);

k0=2*pi/lambda;
d=10*lambda;

if pol==0
  ff=Mu;
else
  ff=Epsilon;
end

gamma=zeros(1,length(Type));

% You may want to change the determination of the square root.
g=length(Type);

gamma=sqrt(Epsilon(Type).*Mu(Type).*k0^2-ones(1,g)*kx^2);
gamma(2:g-1)=gamma(2:g-1).*(1-2*(imag(gamma(2:g-1))<0));
gamma(1)=extsqrt(Epsilon(Type(1))*k0.^2-kx.^2);
gamma(2)=intsqrt(Epsilon(Type(2)).*Mu(Type(2)).*k0^2-kx^2);
gamma(g-1)=intsqrt(Epsilon(Type(g-1)).*Mu(Type(g-1)).*k0^2-kx^2);
gamma(g)=extsqrt(Epsilon(Type(g))*k0.^2-kx.^2);

% First transfer matrix

b1=gamma(1)/ff(Type(1));
b2=gamma(2)/ff(Type(2));
T1=[b1+b2,b2-b1;b2-b1,b1+b2]./(2*b2);

% Scattering matrices for the inside of the waveguide.

for k=1:length(Type)-3
  t=exp(i*gamma(k+1)*hauteur(k+1));    
  T{2*k-1}=[0,t;t,0];    
  b1=gamma(k+1)/ff(Type(k+1));
  b2=gamma(k+2)/ff(Type(k+2));
  T{2*k}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
end
t=exp(i*gamma(length(Type)-1)*hauteur(length(Type)-1));    
T{size(T,2)+1}=[0,t;t,0];

% Last transfer matrix
b1=gamma(length(Type)-1)/ff(Type(length(Type)-1));
b2=gamma(length(Type))/ff(Type(length(Type)));
TN=[b1+b2,b2-b1;b2-b1,b1+b2]./(2*b2);

H{1}=T{size(T,2)};
A{1}=T{1};

  for j=1:size(T,2)-1
	
    A{j+1}= cascade(A{j},T{j+1}); 	
    H{j+1}= cascade(T{size(T,2)-j},H{j});
  end

S=A{length(A)};


% Intermediary coefficients.

  for j=1:size(T,2)-1
	I{j}=[A{j}(2,1)*H{size(T,2)-j}(1,1),H{size(T,2)-j}(1,2);A{j}(2,1),A{j}(2,2)*H{size(T,2)-j}(1,2)]./(1-A{j}(2,2)*H{size(T,2)-j}(1,1));
  end

r=(-TN(1,1)/TN(1,2)-S(2,2))*(T1(1,1)/T1(2,1)-S(1,1))-S(2,1)*S(1,2);

 % Intermediary coefficients outside of the scattering matrices 

I3=[T1(1,1)/T1(2,1),0;1,0];
IN2=[0,1;0,-TN(1,1)/TN(1,2)];


I2=[1/T1(2,1),0;0,0];
IN1=[0,0;0,TN(2,1)-TN(2,2)*TN(1,1)/TN(1,2)];


I1=I2.*exp(i*gamma(1)*hauteur(1));
IN=IN1.*exp(i*gamma(length(Type))*hauteur(length(Type)));

I=[I1,I2,I3,I,IN2,IN1,IN];

% We choose arbitrary one of the unknowns.

B2=1;
A2=T1(1,1)/T1(2,1)*B2;
AN=B2*(T1(1,1)/T1(2,1)-S(1,1))/S(1,2);

%_____________________________________________________________________________
% Computation of the field in the different layers

t=1;
h=0;

E=zeros(sum(ny),1);	%initialisation du veteur

k=1;

for k=1:length(Type)
   for m=1:ny(k)
      h=h+hauteur(k)/ny(k);
	E(t,1)=(I{2*k}(1,1)*B2*exp(i*gamma(k)*(hauteur(k)-h))+I{2*k}(1,2)*AN*exp(i*gamma(k)*(hauteur(k)-h))+I{2*k-1}(2,1)*B2*exp(i*gamma(k)*(h))+I{2*k-1}(2,2)*AN*exp(i*gamma(k)*(h)));

	t=t+1;
	end
    h=0;
end

t=(1:t-1);


% -------------- Vizualisation --------------------------

for i=1:length(hauteur)
y(i,1)=sum(ny(1:i));
end
y(:,2)=y(:,1)+1/(y(length(hauteur),1)*1000);

V1=real(E);

figure(1)
hold on
plot(t',V1(t,1),'linewidth',4),xlabel('Position in the slab thickness'),ylabel('Field');
for i=1:length(hauteur)
m=min(V1(t,1)); M=max(V1(t,1));
plot(y(i,:),[m,M],'r','linewidth',2);
end
hold off


end

end
