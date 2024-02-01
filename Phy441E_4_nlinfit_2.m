rx=1:1:30; % lets create our x vector
ry=1:1:30; % lets create our y vector

[X,Y]=meshgrid(rx,ry); % this creates the mesh grid matrix depending on ry & rx

xc=15;yc=15;A=30;w=5;  % constants for our 3D Gaussian Function
n=[A,w,xc,yc];         % for creating our 3D Gaussian Function

x=X(:);y=Y(:);         % our 3D Gaussian Function accepts a Nx2 matrix 
                       % where (:,1)->x, (:,2)->y      %coluonmized 
M=[x,y];               % Do not forget these are gridded values of X and Y.  
G=gauss2D(n,M);        % now it is time to create the Gaussian Z numbers (aka, color intensities)
g=reshape(G,30,30);    % However G is going to be in Nx2 size , i.e. N=900 here....
figure(13);surf(X,Y,g);hold on; % inorder to make it proper 3D figure, we need to use reshape !!!

%%%%% Creating the DATA ENDS HERE !!!! %%%%%


% Fitting in 3D %%% using nlinfit !!!!!! %
nuve=[2,1,13,18]; % nuve is created !gedit 
memo_options=statset;
memo_options.Display='iter'; % nlinfit uses statset not optimset !!!! 
% regular usage is below:
%[beta2,R2,J2,covB2]=nlinfit(q,A,@gaussfit,n2,'Weights',dA);
% but we need to make some changes... for Nx2 data usage 

Xnew=[x,y]; % our new X matrix which Nx2 matrix
Ynew=G;     % Y values which is still Nx1 matrix
dYnew=sqrt(Ynew); % our new dA errorbar matrix Nx1 


[nfit,R,Jac,covB]=nlinfit(Xnew,Ynew,'gauss2D',nuve,memo_options,'Weights',dYnew);
%[beta,R,J,covB]=nlinfit(Xnew,Ynew,'gauss2D',n,'Weights',sqrt(Ynew));
%nfit = 30.0000    5.0000   15.0000   15.0000

%%% OR USING nlfit3.m 
optim_para.figure_no=111; 
optim_para.mask=[1 1 1 1];
[fv,ev,chi2,fr] = nlfit3(Xnew,Ynew,'gauss2D',nuve,optim_para,dYnew);  


Gfit=gauss2D(nfit,Xnew);                                   % Gfit is created !
Gfit=reshape(Gfit,30,30);                                  % Gfir is put in rx x ry size;        
figure(14);surf(X,Y,Gfit);hold on;

% calculating statistics about the fit %
noffreedom=length(Xnew)-length(nfit);               % or noffreedom=length(Xnew(:))-length(nfit)
ci = nlparci(nfit,R,'jacobian',J);
errobars=(abs(abs(ci(:,1))-abs(ci(:,2))))*0.5;
chi2m=((R'*R)./(sum(dYnew.^2)))*(1/noffreedom) ;     %This will the chi2 values
