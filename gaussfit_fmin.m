function [Afit]=gaussfit_fmin(n)
% simple Gaussian definition for a Gaussian fit. 
% beta0=[Amplitude,Centre,Width,AConstantBackground];
% this is for fminsearch example for searching the Peak CENTER parameter !
% at this moment it is for lc4O8_d0025_a , scan 29 !!! 

Xc=n(1);

Amp=1500;
width=0.01;
backgnd=1500;

load q_mat.mat;

%Afit=(n(1)*exp(-1*(X-n(2)).^2/(n(3)^2))+n(4));
Amodel=(Amp*exp(-0.5*(q-Xc).^2/(width^2))+backgnd);
Afit=max(A(:)-Amodel(:));
     %in order to search the peak parameter best value, we can use the max of the value above!
