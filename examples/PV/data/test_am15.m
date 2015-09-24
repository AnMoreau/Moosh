
more off

N=4000-280+1;
x=linspace(375,750,N);
y=x;
for j=1:N
	y(j)=am1_5(x(j));
end

	plot(x,y)
	power=sum(y)*(max(x)-min(x))/N
	power=trapz(x,y)

clear
  N=500
  x=linspace(375,750,N);
y=x;
for j=1:N
	y(j)=solar(x(j));
end
figure(2)
  plot(x,y)
  current=sum(y)*(max(x)-min(x))/N
  current=trapz(x,y)
