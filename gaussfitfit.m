function [Afit]=gaussfitfit(x,a,b,c,d)
% simple Gaussian definition for a Gaussian fit. 
% beta0=[Amplitude,Centre,Width,AConstantBackground];

Amp=abs(a);
Xc=(b);
width=abs(c);
backgnd=abs(d);


%Afit=(n(1)*exp(-1*(x-n(2)).^2/(n(3)^2))+n(4));
Afit=(Amp*exp(-0.5*(x-Xc).^2/(width^2))+backgnd);
