function [Afit]=gaussfitconv(n,X,Ar)
% written at ITU on Dec 2023 by Mmt
Amp=abs(n(1));
Xc=n(2);
width=abs(n(3));
backgnd=abs(n(4));
Agauss=(Amp*exp(-0.5*(X-Xc).^2/(width^2))+backgnd);

Arn=Ar./max(Ar);                               %noryap=max(Ar)/max(Agauss);
Acon=conv(Agauss,Arn,'same');                  %Afit=Acon./noryap;

Afit=Acon;
