% This is the very much proper fit for a Scattering Experiment Peak involving 
% Resolution Scan and its convalution 
% Dec 2023 at ITU , by MMT

[verim]=specokuma('lc4O8_d0025_a',29);
[verimr]=specokuma('lc4O8_d0025_a',28);
qr=verimr(:,1);Ar=verimr(:,end);dAr=sqrt(verimr(:,end));
q=verim(:,1);A=verim(:,end);dA=sqrt(verim(:,end));
figure; hold on;

% the length pf data is : 51 while the length of resolution is 81 .
% before we use resolution scan , we need it to be normalized ! 
% or we normalized it in fitting routine. 

errorbar(q,A,dA,'ok');
%Ac=convn(A,Ar,'same');  % since we use 'same' we do not need padding ! in theory !
%noryap=max(Ac)/max(A);
%errorbar(q,Ac./noryap,dA,'ob');
n2=[3000 0.22 0.01 100]; 
mask=[1 1 1 1]; 
[nfit,chisquare,errors,fitresult] = nlfit00_class(q,A,'gaussfitconv',n,mask,dA,Ar);

figure; hold on;
errorbar(q,A,dA,'ok');
l1=plot(q,fitresult.yfit,'g-');
% but if you see the figure does not look good in the left and right ends.
%meaning PADDING is needed for a convolution routine. 

Aleft=ones(25,1)*A(1);
Aright=ones(25,1)*A(end); 
Am=[Aleft;A;Aright];
for ilk=1:25
 qleft(ilk)=q(1)-ilk*0.002;
 end;
qlft=qleft(end:-1:1);

for ilk=1:25 
 qright(ilk)=q(end)+ilk*0.002;
 end;
 
qm=[qlft';q;qright'];
dAm=sqrt(Am);

[nfit,chisquare,errors,fitresult] = nlfit00_class(qm,Am,'gaussfitconv',n,mask,dAm,Ar);

yfit=fitresult.yfit(26:76);
l2=plot(q,yfit,'r-');

nfit =
       193.71
      0.22456
    0.0078251
        157.7

[nfit1,chisquare1,errors1,fitresult1] = nlfit00_class(q,A,'gaussfit',n,mask,dA);
nfit1 =
       1416.6
      0.22691
     0.016131
         1664
         
%So the corrected width is 0.0078 not 0.0161 !!!!! 
 
