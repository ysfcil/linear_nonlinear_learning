function [Amin]=gaussfit_lsqnonlin_v1(n,q,A)
% simple Gaussian definition for a Gaussian fit. 
% beta0=[Amplitude,Centre,Width,AConstantBackground];
% This is SPECIFIALLY for Phy508E_6_lsqnonlin_Stuff.m 
%where we seek the min of the (Ydata-YModel)min for the best fit parameters. 

Amp=abs(n(1));
Xc=n(2);
width=abs(n(3));
backgnd=abs(n(4));

%load q_mat.mat; %since we cannot get the A and q values, into the gaussfit_lsqnonlin function
                %we can only cheat the requirements by loading the data externally.

Afit=(n(1)*exp(-1*(q-n(2)).^2/(n(3)^2))+n(4));
%Afit=(Amp*exp(-1*(q-Xc).^2/(width^2))+backgnd);
Amin=Afit-A;    % A & q values are read from q_mat.mat file for lc4O8_d0025_a , scan 29 !!! 
