% Our Data in this subject is an X-ray Diffraction Area Data collected in 2012 
% Argon National Lab. Chicago. The data is a surface of a quasi-crystal approximate sample
% the resolution for MAR-Area detector is : 2300 pixel in both X and Y direction. 
%file name: crru_101_279.raw2300

res=2300;
okunacak='crru_101_279.raw2300';
fid=fopen(okunacak,'r');
I=fread(fid,[res,res],'uint16'); % fread is used to collect data from the data file.
fclose(fid); % the file is written for un-sign int 16 type. 

A=rot90(I); % the data is needed to be turned by 90 deg to put into correct direction.
figure(1);surf(A);view(2);colormap(bone);shading flat;axis square; caxis([0 50]);colorbar;
axis('square');axis([0 2300 0 2300]);

%now we see our data. The focus goes to the left white spot on the left side of centre.
%so we need to focus on this spot.

Afo=A([1000:1200],[600:800]); % The matrix is designed as : A[y,x] infact;
figure(2);surf(Afo);view(2);colormap(bone);shading flat;axis square;caxis([0 50]);
axis('square');


% now lets do this focusing a little bit better.
Afo=A([1100:1200],[700:800]); % here we can also use view([xmin xmax ymin ymax]) but this only
                              % focus the view , not focus the data.
figure(3);surf(Afo);view(2);colormap(bone);shading flat;axis square;caxis([0 50]);
axis('square');axis([0 100 0 100]);colorbar; 
% we check how this looks in full c-axis scale as well 
caxis('auto'); 
caxis([0 100]);

%now we need to create of x & y axis in the same amount with our focused data.
size(Afo); %-> 101x101 

x=1:1:101;
y=1:1:101;
[xgrd,ygrd]=meshgrid(x,y);
Xn=[xgrd(:),ygrd(:)]; % we need to send a X matrix and corresponding Y values to be fit.
Aduz=Afo(:); % Therefore, we need to linearized the Y values and make X as (N:,1)->x
             %                                                        Y as  (N:,2)->y
Aduz      10201x1              81608  double              
Xn        10201x2             163216  double

%Our 3D Gaussian function is listed under the name of : fitAlan.m
%%bg=n(1);Amp=n(2);xc=n(3);wx=abs(n(4));yc=n(5);wy=abs(n(6));
n=[10 1e4 40 5 50 5 ]; % center is roughly at (40,50) !!!
[Abas]=fitAlan(n,Xn); % check how our starting looks like.

AbasMat=reshape(Abas,size(Afo)); %now we put our calculated Z values for (x,y) in the same
                                 % size of the original data. 
figure(5);surf(AbasMat);view(2);colormap(bone);shading flat;axis square;caxis('auto');
axis('square');axis([0 100 0 100]);colorbar;

%Since we created a 3D matrix using the 3D gaussian definition in fitAlan.m
% it is just the right time to start fitting.
param=statset('Display','iter'); 
%options=optimset('Display','iter','TolFun',1*10^(-9)); % this would be enough for nlfit !!!
[pfit,r,J,con,mse]=nlinfit(Xn,Aduz,'fitAlan',n,param);

pfit =
        23.67        68660       38.187        1.515       49.074       1.2531

% now first 1) create our fit picture, 2) calculate ChiSquare and 3) errorbars.
[Afit]=fitAlan(pfit,Xn); % this puts Afit into a (Nx1) matrix so we need it reshape
AfitMat=reshape(Afit,size(Afo)); %using the original data matrix Afo : Focused data
figure(6);surf(AfitMat);view(2);colormap(bone);shading flat;axis square;caxis([0 100]);
axis('square');axis([0 100 0 100]);colorbar; 

noffreedom=length(Xn(:))-length(pfit); %or noffreedom=length(Xn(:,1))-length(nfit);
ci = nlparci(pfit,r,'jacobian',J);

ci =
       21.419       25.921
        68476        68844
       38.184        38.19
       1.5109       1.5191
       49.072       49.077
       1.2498       1.2564


hatabars=(abs(abs(ci(:,1))-abs(ci(:,2))))*0.5 ; % errorbar lar burda !!!  
chi2m=((r'*r)./(sum(dAduz.^2)))*(1/noffreedom); 



% now if we want to fit with errorbars included !! Weighted in other words.
dAduz=sqrt(Aduz);
[pfit1,r1,J1,con1,mse1]=nlinfit(Xn,Aduz,'fitAlan',n,options,'Weights',dAduz);


 
