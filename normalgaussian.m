function y=normalgaussian(x)
%normal gaussian : A is here Area and omega:2sigma=0.849FWHM
% in other words, this is an AREA GAUSSIAN !!! 

A=100;
xc=0;
omega=5;
f1=A/(sqrt(pi/2)*omega);
f2=exp(-2*((x-xc)/omega).^2);
y=f1.*f2;

