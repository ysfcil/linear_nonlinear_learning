% this is for the nonlinear teaching material created on Nov 2023 at ITU by Mmt.

load census;
options = fitoptions;
options.Normalize = 'on';  % options=fitoptions('Normalize','on');
figure(1); plot(cdate,pop,'ob'); hold on;

f1=fit(cdate,pop,'poly3',options); % fitting with 3rd order polinomial
f2=fit(cdate,pop,'exp1',options);
f3=fit(cdate,pop,'cubicspline',options)
 
% example : >>f1
              f1 = 
              Linear model Poly3: f1(x) = p1*x^3 + p2*x^2 + p3*x + p4
                                  where x is normalized by mean 1890 and std 62.05
                                  Coefficients (with 95% confidence bounds):
                                  p1 =       0.921  (-0.9743, 2.816)
                                  p2 =       25.18  (23.57, 26.79)
                                  p3 =       73.86  (70.33, 77.39)
                                  p4 =       61.74  (59.69, 63.8)
                                  
plot(f1);plot(f2);plot(f3);  % directly plots the fits !!!

[f,gof,out] = fit(cdate,pop,'SmoothingSpline');  % fitting with SmoothingSpline Polynomial !
              % f: is the fit >>plot(f,'-g');
              % gof= 
              struct with fields:
              sse: 11.391
              rsquare: 0.99991
              dfe: 7.751
              adjrsquare: 0.99976
              rmse: 1.2123
              %out= 
              struct with fields:
              numobs: 21
              numparam: 13.249
              residuals: [21x1 double]
              Jacobian: []
              exitflag: 1
              p: 0.0089197
% So we can continue to calculate chi2 as;
r=out.residuals;
(r'*r)/out.numparam
ans =
      0.85976

options=fitoptions('Method','SmoothingSpline',...
                     'SmoothingParam',0.0098);
[f,gof,out] = fit(cdate,pop,'SmoothingSpline',options);
plot(f,'-g');
    
    
%%%%%%%%%%%%%%%%%%%%% Gaussian DATA fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
a1 = 1; b1 = -1; c1 = 0.05;
a2 = 1; b2 = 1; c2 = 50;
x = (-10:0.02:10)';
gdata = a1*exp(-((x-b1)/c1).^2) + ...
        a2*exp(-((x-b2)/c2).^2) + ...
        0.1*(rand(size(x))-.5);
plot(x,gdata); hold on;

gfit = fit(x,gdata,'gauss2');
 
% This will not fot the data !!!
options = fitoptions('gauss2', 'Lower', [0 -Inf 0 0 -Inf 0]);
             % or options = fitoptions('gauss2');
             %    options.Lower = [0 -Inf 0 0 -Inf 0];
gfit = fit(x,gdata,'gauss2',options); 

options = 
  nlsqoptions with properties:

       StartPoint: []
            Lower: [0 -Inf 0 0 -Inf 0]
            Upper: []
        Algorithm: 'Trust-Region'
    DiffMinChange: 1.0000e-08
    DiffMaxChange: 0.1000
          Display: 'Notify'
      MaxFunEvals: 600
          MaxIter: 400
           TolFun: 1.0000e-06
             TolX: 1.0000e-06
           Robust: 'Off'
        Normalize: 'off'
          Exclude: []
          Weights: []
           Method: 'NonlinearLeastSquares'

%%%%%%%%%%% WRITTING YOUR FIT LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lft = fittype({'x','sin(x)','1'});
% lft = 
%     Linear model:
%     lft(a,b,c,x) = a*x + b*sin(x) + c
fo = fitoptions(lft)
fo = 
  llsqoptions with properties:

        Lower: []
        Upper: []
       Robust: 'Off'
    Normalize: 'off'
      Exclude: []
      Weights: []
       Method: 'LinearLeastSquares'

fo.Normalize = 'on';

%%% The list of the Model Functions ::::
 https://www.mathworks.com/help/curvefit/list-of-library-models-for-curve-and-surface-fitting.html#btbcvnl
 
%%%% Fitting our Data !!!! %%%%%%%%%%%%
% Now please do not forget to go to ./non-linear-A directory !!! 
[verim]=specokuma('lc4O8_d0025_a',29); 
q=verim(:,1);A=verim(:,end);dA=sqrt(verim(:,end)); 
figure(3); hold on;
errorbar(q,A,dA,'ob'); 
% bacause we cannot simply >> lmodel=fittype({'x','gauss1','1'}); % We cannot !!!! 
% edit our below gaussfit.m function for gaussfitfit.m as below:
                         %function [Afit]=gaussfitfit(x,a,b,c,d)
                         % simple Gaussian definition for a Gaussian fit using fit command. 
                         % beta0=[Amplitude,Centre,Width,AConstantBackground];

                         %Amp=abs(a);
                         %Xc=(b);
                         %width=abs(c);
                         %backgnd=abs(d);

                         %Afit=(Amp*exp(-0.5*(x-Xc).^2/(width^2))+backgnd);
lmodel = fittype( 'gaussfitfit(x,a,b,c,d)' );
f = fit( q, A, lmodel, 'StartPoint', [3000 0.22 0.01 100]); 
f = 
     General model:
     f(x) = gaussfitfit(x,a,b,c,d)
     Coefficients (with 95% confidence bounds):
       a =        1422  (1314, 1531)
       b =      0.2263  (0.2251, 0.2274)
       c =     0.01527  (0.01363, 0.01692)
       d =        1698  (1617, 1779)
plot(f,'-k');
  
%%% adding Weigths and more options %%%%%% 
opsiyon=fitoptions(lmodel)              
opsiyon =
        Normalize: 'off'
          Exclude: []
          Weights: []
           Method: 'NonlinearLeastSquares'
           Robust: 'Off'
       StartPoint: [1x0 double]
            Lower: [1x0 double]
            Upper: [1x0 double]
        Algorithm: 'Trust-Region'
    DiffMinChange: 1e-08
    DiffMaxChange: 0.1
          Display: 'Notify'
      MaxFunEvals: 600
          MaxIter: 400
           TolFun: 1e-06
             TolX: 1e-06

opsiyon.Weights=dA;
opsiyon.StartPoint=[3000 0.22 0.01 100];
opsiyon =
        Normalize: 'off'
          Exclude: []
          Weights: [40.608 40.62 40.012 39.421 40.915 40.05 40.237 40.645 40.534 40.485 40.299 ... ]
           Method: 'NonlinearLeastSquares'
           Robust: 'Off'
       StartPoint: [3000 0.22 0.01 100]
            Lower: [1x0 double]
            Upper: [1x0 double]
        Algorithm: 'Trust-Region'
    DiffMinChange: 1e-08
    DiffMaxChange: 0.1
          Display: 'Notify'
      MaxFunEvals: 600
          MaxIter: 400
           TolFun: 1e-06
             TolX: 1e-06

f1 = fit( q, A, lmodel,opsiyon);
     f1 = 
     General model:
     f1(x) = gaussfitfit(x,a,b,c,d)
     Coefficients (with 95% confidence bounds):
       a =        1420  (1313, 1528)
       b =      0.2259  (0.2249, 0.227)
       c =     0.01486  (0.01325, 0.01646)
       d =        1717  (1634, 1800)
% this is abit different than f (where we do not have dA errorbars in the fitting )! 

  
