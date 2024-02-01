function y=gaussprimitive(x)
sig=2.5;
A=100;
xc=0;
f=A*exp(-(x-xc).^2/sig^2);
y=f;


