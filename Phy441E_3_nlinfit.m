% this version is updated on March 2023 , at ITU...

[verim]=specokuma('lc4O8_d0025_a',29); 
q=verim(:,1);A=verim(:,end);dA=sqrt(verim(:,end)); 
figure(3); hold on;
errorbar(q,A,dA,'ob');
nuve=[3000 0.22 0.01 100];                                   
[nfit1,R1,J1,covB1]=nlinfit(q,A,@gaussfit,nuve,'Weights',dA);
%[nfit,R1,J1,covB1]=nlinfit(q,A,'gaussfit',nuve,'Weights',dA);
%statset is my option table variable whoch we need to tweak some parts of it to see the details
opsiyon=statset('Display','iter');
%[nfit,R2,J2,covB2]=nlinfit(q,A,'gaussfit',nuve,opsiyon);
[nfit2,R2,J2,covB2]=nlinfit(q,A,'gaussfit',nuve,opsiyon,'Weights',dA);
%tic;sum(R2.^2);toc                  
%Elapsed time is 0.000340 seconds.
%tic;(R2'*R2);toc  
%Elapsed time is 0.000167 seconds.
% So because the second way is faster than the first way, so we folow the second way . 
% This is the way !!!!  

Afit=gaussfit(nfit,q);
plot(q,Afit,'-r');

noffreedom=length(q)-length(nfit); % This is our new num of freedom.
ci1 = nlparci(nfit1,R1,'jacobian',J1); %ci = nlparci(beta,R,'jacobian',J); % this will give the up down confidence regime

errobars1=(abs(abs(ci1(:,1))-abs(ci1(:,2))))*1; % the amplitude is a choice ;
errobars1=(abs(abs(ci1(:,1))-abs(ci1(:,2))))*0.5; % this is much more correct one !  

ci2 = nlparci(nfit2,R2,'jacobian',J2); %ci = nlparci(beta,R,'jacobian',J); % this will give the up down confidence regime
errobars2=(abs(abs(ci2(:,1))-abs(ci2(:,2))))*1; % the amplitude is a choice 

%chi2m=((R2'*R2)./(sum(dA.^2)))*(1/noffreedom) ; %This will the chi2 values
chi2_1=((R2'*R2)./(sum(dA.^2)))*(1/noffreedom)
% chi2_1=((R1'*R1)./(sum(dA.^2)))*(1/noffreedom) -> this is a bit wrong

chi2_1=(R1./dA)'*(R1./dA)/noffreedom  % 8.28
chi2_1=(R2./dA)'*(R2./dA)/noffreedom  % 397.4 
 
%%%(error_bar2=full((1/1).*sqrt(((Cii2)).*rnorm2/(length(X)))*tinv(0.84135,noffreedom2))


