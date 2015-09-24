% 
%   [z,msg]=descent(lambda,kx_start,step,precision,@f)
% 
% Makes a steepest descent on "f" in order to find its closest zero.
% Stops when the value of "f" is smaller than "precision"
% The initial value of the step is "step" - increase it if the descent
% is too slow or if it doesn't reach a solution in less than 1000 steps
% (in which case it returns 1000 in "msg", meaning no solution was found).

function [z,msg]=descent(lambda,z0,step,precision,fct)

% Initialization
z=z0;
a=fct(z,lambda);
j=0;
msg=0;
% We take arbitrarily a step to compute the gradient that is 100 times
% smaller than the descent step.
dz=0.01*step;

while ((a>precision)&&(j<1000))
% Computation of the gradient
grad=(fct(z+dz,lambda)-fct(z-dz,lambda))/(2*dz)+i*(fct(z+i*dz,lambda)-fct(z-i*dz,lambda))/(2*dz);
  if (abs(grad)~=0)
    grad=grad/abs(grad);
  else
    % Hum, why would we be exactly zero ? That's weird...
    % Let's make it an error
    msg=1000;
    break
  end

  z_new=z-step*grad;  
  b=fct(z_new,lambda);
  if (b>a)
    step=step/2.;
    dz=dz/2.;
  else
    a=b;
    z=z_new;
  end

  j=j+1;

end

if (msg~=1000)
  msg=j;
end

end
