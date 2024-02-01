function [Afit]=gaussfit(n,X)
% simple Gaussian definition for a Gaussian fit. 
% beta0=[Amplitude,Centre,Width,AConstantBackground];

Amp=abs(n(1));
Xc=n(2);
width=abs(n(3));
backgnd=abs(n(4));


%Afit=(n(1)*exp(-1*(X-n(2)).^2/(n(3)^2))+n(4));
Afit=(Amp*exp(-0.5*(X-Xc).^2/(width^2))+backgnd);
