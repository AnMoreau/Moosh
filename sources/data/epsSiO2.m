%ordinary index of SiO2 from Ghosh 1999 a-quartz n(o).

function n2=epsSiO2(lambda)
lambda=lambda*1e-3;
n2=1.28604141+(1.07044083*lambda.^2)/(lambda.^2-1.00585997e-2)+(1.10202242*lambda.^2)/(lambda.^2-100);
end
