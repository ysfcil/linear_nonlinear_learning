% Let's create our data (linear - like) and the figure for it.
x=0:1:7;
y=[0,20,60,68,77,110,100,130];
dy=sqrt(y);

figure(1); hold on;
errorbar(x,y,dy,'bo');

[nfit,others]=polyfit(x,y,1); % background and slope
yfit=polyval(nfit,x);

plot(x,yfit,'-r');  % this is the fit 

yman=nfit(2)+nfit(1).*x; % nfit(2) is the constant and nfit(1) is the slope.
[yfit1,delta] = polyval(nfit,x,others); % others is struct .... We have estimated erros 
%%%%%% Solving with mldivide !!!! %%%%%%%
slope=x'\y';                      % we need y=mx for each point so take tranpose !
yslope=slope*x;                   % now create the fit line of y=mx;
plot(x,yslope,'-g');  % this fits using mldivide command a simple y=mx ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Fitting with y=bx+c with mldivide !! %%%%%%%
xnew=[ones(length(x),1) x']; % we create y=1+mx for each data point ;
slope2=xnew\y'               % this fits y=bx+c; 
    slope2 =
       9.4167
       17.488      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% ---- The Results Of the Fit Analysis ----
others.normr=>28.37;
others.df=>6; %(because N=8, fitting parameters are in the number of m=2; (back and slope)
chi2=(others.normr/others.df); 

% because; others.norm is actually;
%  equals to; 
% N=8 here , and m=2 : N-m : dimension of freedom, freedom number
  residu=sum(abs(y-yfit).^2);
->804.87
  norm_residu=sqrt(residu);
->28.37

% and the errorbar on calculated fitting parametres are:

R=others.R; % Residue matrix
cov=(R*R')*others.normr.^2/others.df ; % covariance matrix
delta=1./sqrt(diag(cov))*tinv(0.84413,others.df); % here (1+.68)/2=0.84 is used.

% so the errorbar(1)=delta(1); errorbar(2)=delta(2);


