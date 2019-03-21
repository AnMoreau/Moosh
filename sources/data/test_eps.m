lam=linspace(200,4000,200);

for j=1:length(lam)

epsi(j)=epsNibb(lam(j));

endfor

plot(-log(lam),log(-real(epsi)),'g','linewidth',2,-log(lam),log(imag(epsi)),'linewidth',2)