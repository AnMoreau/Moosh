function S=couche(b,h,Kl)
t=exp(+i*b*h);
l=exp(-1*Kl*h);
S=[0 0 t 0; 0 0 0 l;t 0 0 0;0 l 0 0];
end



    
