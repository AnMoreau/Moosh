%
% S=cascade(A,B)
% This function takes two 2x2 matrices A and B, that are assumed
% to be scattering matrices and combines them assuming A is the "upper"
% one, and B the "lower" one, physically. The result is a 2x2 
% scattering matrix.
%

function S=cascade(A,B)
t=1/(1-B(1,1)*A(2,2));
S=[A(1,1)+A(1,2)*B(1,1)*A(2,1)*t,A(1,2)*B(1,2)*t;B(2,1)*A(2,1)*t,B(2,2)+A(2,2)*B(1,2)*B(2,1)*t];
end