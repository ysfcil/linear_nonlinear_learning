% Written by Mmt on July 2022 at Emirgan, Ist.
% The example uses the fminsearch command !
% fminsearch cammand uses only a number not a matrix. So we can only fit a parameter
% best for our needs. Afit-> a number !!!!
% uses gaussfit_fmin.m at (Ydata-YModel)min
% edit March 2023, ITU for gauss_fmin_v1 better and more proffesional usage 

[verim]=specokuma('lc4O8_d0025_a',29); 
q=verim(:,1);A=verim(:,end);dA=sqrt(verim(:,end)); 
figure(3); hold on;
errorbar(q,A,dA,'ob');

%nuve=[3000 0.22 0.01 100]; 
nuve=0.22;

nfit=fminsearch('gaussfit_fmin',nuve);

 %results: nfit = 0.22655 !!! 
   %this is fairly robust result since nuve=0.18 and nuve nuve=0.25 gives the same result !
   
options = optimset('Display','iter');
nfit=fminsearch('gaussfit_fmin',nuve,options);
   %gives the iterations !

%%options = optimset('Display','iter','PlotFcns',@optimplotfval);

nfit=fminsearch('gaussfit_fmin_v1',nuve,options,q,A); % now you do not need unncessary
                                                      % qmat.mat to be read in the model. 
