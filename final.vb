
%Yusuf Cil 090200118
function [beta,chisquare,errors,fitresult] = nlfit00_class(X,y,model,beta0,iterprint, iterplot,varargin)
%This is the non-linear fit code using Marquardt-Levenberg Tech. 
%This is my final project for FIZ 441E CRN:11080

%now if any optional parameters are specified, set a and sigma to the 
%appropriate values



linetype = ["r-","r-","r-","r-","r-", "r--","r:","r-.","g-","g--","g:","g-.",
 "b-","b--","b:","b-.", "c-", "c--","c:","c-.","m-","m--","m:","m-.", "y-","y--","y:","y-.", "k-",
 "k--","k:","k-.","w-","w--","w:","w"];
%first ones have the samestyling because they will be far from the data points
switch nargin  
case 7
   a=deal(varargin{1});
case 8
   [a,sig]=deal(varargin{:});
   if isempty(a)
      a=ones(size(beta0)); %default to fitting all parameters
   end
case 9
    [a,sig,fun_struct]=deal(varargin{:}); %any additional parameters are passed on to model function 
    if isempty(a)                   %as a structure
        a=ones(size(beta0));
        sig=sqrt(y(:));
    end
end

 
n = length(y);  %define number of pts
if min(size(y)) ~= 1
   error('Requires a vector second input argument.');
end
X = X(:);  %columnize vectors
y = y(:);
beta0 = beta0(:);  
sig=sig(:);
sig(sig==0)=inf;  %guard against division by zero

A=find(a(:));  %A defines which parameters are being varied. eg A=[1 2 4];
p=length(A);   %define number of parameters being fitted.

J = zeros(n,p);    %initialize Jacobian
beta = beta0;      %initialize parameters
betanew=beta0 * 1.01;   %slightly vary the parameters.  A serves to ...
                          %pick out only those parameters being varied.   
betaLM=beta0;

iter = 0;
tolerance = 1e-4; %fractional tolerance of variations in parameters and fractional tolerance on chi-square
sse = 1;     %sse is Chi-square
sseold = 1;
lambda=0.01; %Marquardt - Levenberg Step size to begin !
maxiter = 39;%maximum allowed iterations
while (~isempty(find((abs((betanew-beta)./(beta+sqrt(eps))) > tolerance), 1)) | ...
      ((sseold-sse)/(sse+sqrt(eps)) > tolerance)) & ...
      (iter < maxiter)
   %while either the fractional change in any parameter beta or chi-square is 
   %greater than specified tolerances AND the number of iterations 
   %is less than max
iter = iter + 1;


if nargin > 8
   yfit = feval(model,beta,X,fun_struct);   %evaluate the model
else
    yfit=feval(model,beta,X);
end
   r = (y - yfit)./sig;          %calculate residual
   sseold = r'*r;                %calculate chi-square

if iterprint     %I put the printing option as argument to the function so the user can
                   %choose to print each iteration or not. It still prints the last one
   fprintf('Iteration #%d \t Chi-Square %f \n',iter,sse/(n-p)); %we display some info for each iteration
end
   
   for k = 1:p
      %calculate J=dy/da
      if nargin > 8
          J(:,k)=nlfit_deriv(model,beta,X,yfit,A(k),nargin,fun_struct)./sig;
      else
          J(:,k)=nlfit_deriv(model,beta,X,yfit,A(k),nargin)./sig;
      end
      %A(k) cleverly passes wrt which parameter in beta should the derivative
      %be calculated as k cycles through the total number of varied parameters. 
      %This step calculates the derivative of the model, yfit, wrt to small changes
      %in the parameters, beta. It calculates it at each point, returning #pts 
      %derivatives.
    
      
   end


   % Levenberg-Marquardt type adjustment 
   % Gauss-Newton step -> J\r
   % LM step -> inv(J'*J+constant*eye(p))*J'*r

   %Gradient search

   JJ=J'*J;  %calculate J squared which is equal to alpha_ll
   Jr=J'*r;  %calculate dChi-square/dbeta(n) = beta						
   stepLM = (JJ-diag(diag(JJ))+JJ.*(eye(p)*(1+lambda)))\Jr; 
   %first step replaces the diagonal elements of JJ, which are equal to alpha_ll
   %with (1+lambda)*alpha_ll
   %Dividing by Jr, which is equal to beta, then computes the next step of the parameters
   %according to LM algorithim. c.f. 15.5.14 in Numerical Recipes in C.

   sseLM = sseold*2; %give sseLM some value so we enter the while loop
   iter1 = 0;        %this saves us from repating lines of code outside
                     %the loop and doesn't change the result

   while sseLM > sseold & iter1 < 12   %keep stepping until we reduce chi-square
      stepLM = stepLM/sqrt(10);
      betaLM(A) = beta(A) + stepLM;
      if nargin > 8
          yfitnew = feval(model,betaLM,X,fun_struct);
      else
          yfitnew = feval(model,betaLM,X);
      end
      rnew = (y - yfitnew)./sig;
      sseLM = rnew'*rnew;
      iter1 = iter1 + 1;
   end

    if iter1 < 12
      lambda=lambda/2;
      betanew=betaLM;
      sse=sseLM;
	else 
      lambda=lambda*10.; %change the lambda until we get the point of lowering the chi square
    end
    beta = betanew;

    if iterplot %Also added a iterplot option to plot the model from each iteration if needed.
    figure(10);hold on;
    plot(X,y,'ob');
    plot(X,yfit, linetype(iter));
    end
end  

if iter == maxiter
   disp('NLINFIT did NOT converge. Returning results from last iteration.');
end
if iterprint == 0
    fprintf("Chi-square is %f after iteration number #%d  \n", sse/(n-p), iter)
end


%main fitting routine is over. now prepare output values
chisquare=sse/(n-p);   %normalize chi-square by degrees of freedom
fitresult.iterations = iter;
fitresult.phi = sse;
fitresult.lambda = lambda;  %prepare an optional structure fitresult
fitresult.NumPts = n;
fitresult.NumParam = p;
fitresult.yfit = yfit;
fitresult.residual = r;
fitresult.sigma = sig;
fitresult.J=J;
errors=nlparci(beta(A),r,J,A); %calculate errors in fit only for parameters being fit






function y=nlfit_deriv(model,beta,X,Y,n,nargin_nlfit,fun_struct)
%nlfit_deriv(beta,x,n): Take the first derivative of y with respect to beta(n) at
%			points X. This is a subroutine of nlfit.
%        
%        nlfit_deriv only calculates dy/dbeta at each point. Later in the main 
%        routine they will all have to be summed and weighted together.

delta = zeros(size(beta));  %initialize delta to zero
delta(n) = sqrt(eps)*beta(n) + eps; %instead of checking zeroes and adding an eps just added an eps
                                    %doesn't change the results and requires one less if in a looped section


if nargin_nlfit > 8
    y1 = feval(model,beta+delta,X,fun_struct);
else
    y1 = feval(model,beta+delta,X);
end
y = (y1 - Y)/delta(n);  %calculate dy/dbeta

while(sum(y)==0 && delta(n)<0.01*beta)
        delta(n)=delta(n)*10.;
        if nargin_nlfit > 7
            y1=feval(model,beta+delta,X,fun_struct);
        else
            y1=feval(model,beta+delta,X);
        end
        y=(y1-Y)/(sqrt(eps)*beta(n));
end

%end of nlfit_deriv




function delta = nlparci(x,f,J,A)
%code with slight modifications from function nlparci.m
%returns one gaussian standard deviation (confidence interval
%equals 68%) estimate of errors on fitted parameters.

%NLPARCI Confidence intervals on parameters of nonlinear models.
%   CI = NLPARCI(X,F,J) returns the 95% confidence interval CI
%   on the nonlinear least squares parameter estimate X, given the 
%   residual sum of squares, F, and the Jacobian matrix ,J, at the solution.
%
%   The confidence interval calculation is valid for systems where 
%   the number of rows of J exceeds the length of X. 
%
%   NLPARCI uses the outputs of NLINFIT for its inputs.
%   Example:
%      [x,f,J]=nlinfit(input,output,model,xinit);
%      ci = nlparci(x,f,J);
%
%   See also NLINFIT.
%

%   Bradley Jones 1-28-94
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.7 $  $Date: 1998/07/10 14:45:54 $

%initialization
if nargin < 3
   error('Requires three inputs.');
end

f = f(:);
[m,n] = size(J);
if m <= n
   error('The number of observations must exceed the number of parameters.');
end

if length(x) ~= n
   error('The length of x must equal the number of columns in J.')
end

% approximation when a column is zero vector
temp = find(max(abs(J)) == 0);
if ~isempty(temp)
   J(temp,:) = J(temp,:) + sqrt(eps);
end

%calculate covariance
[~, R] = qr(J,0);
Rinv = R\eye(size(R));
diag_info = sum((Rinv.*Rinv),2);

v = m-n;
rmse = sqrt(sum(f.*f)/v);

% calculate one standard deviation
%note that the number .8413 is .5(1+.68), the proper
%formula for calculating a 68% confidence integral
%using the t-distribution.

delta(A) = sqrt(diag_info) .* rmse*tinv(0.8413,v); 

%--end of nlparci.m---