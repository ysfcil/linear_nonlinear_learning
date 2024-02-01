First, we will create a numerical ditribution which is a simle Gauss one. This will be our data points. We also need some errorbars to paly with. Therefore, 

X=1:1:100; % This the x numbers , independent variables stricly speaking. 
A=1000;Xo=50;W=5; % As you may guess, Amplitude is 10; the central position , peak position aka
                % is at 50 and the width is 5. Here we assume bck=0 : background is zero.
Y=A.*exp(-1*(X-Xo).^2/W.^2); % a simple Gaussian , but again strickly speaking this is not
                               % a good way of creating a Gauss.. We will come this point 
                               % later on... Let's go deep into the usage of the code at the 
                               % moment
dY=rand(100,1).*2; % dY is our errorbar matrix. So far we have X , Y, and dY which are all
                               % 1x100 matrices.

figure(1); errorbar(X,Y,dY,'bo'); % first let's see the graph of it.

%Then let's create what we need now in order to run the lsqcurvefit code.
% [nfit,rnorm,r,ef,oput,L,jac]=lsqcurvefit(@gaussfit,bet0,X,Y,lowerb,upperb,options);
% We have X and Y, we do not have MODEL yet. And beta0 is   
%the starting point of the parameters which are going to be fit. What are they ? 
% They are A, Xo, and W at the moment. Remember we do not have backgorund as defined. Or in %
%other words, it is zero for now. 

nuve=[10 50 1 1]; % For God's sake, let's cheat now by staring the parameters very close to %
                  %the real ones. 
%Now let's write the code below: and Save it as LSGauss.m , OK?

function [Afit]=gaussfit(n,X)
% simple Gaussian definition for a Gaussian fit. 
% beta0=[Amplitude,Centre,Width,AConstantBackground];

Amp=abs(n(1));
Xc=n(2);
width=abs(n(3));
backgnd=(n(4));


%Afit=(n(1)*exp(-1*(X-n(2)).^2/(n(3)^2))+n(4));
Afit=(Amp*exp(-1*(X-Xc).^2/(width^2))+backgnd);

Then with the code line below; try to fit it...
>>[nfit1,rnorm,r,ef,oput,L,jac]=lsqcurvefit(@gaussfit,nuve,X,Y,lowerb,upperb,options);

However we still do not have the lowerb and upperb and options...
>> options=optimset('Display','iter-whos
detailed','MaxFunEvals',10000,'MaxIter',5000,'TolFun',1*10^(-4));

lowerb= [] ;
upperb=[];
% to get the details of the optimset , just type 
>> help optimset

When you type the fit command, you will get the fit values, itterations etc... displayed on the screen.
Now we have nfit,rnorm,r,ef,oput,L,jac
nfit : the fitted parameters : It is a vector.
rnorm: squared residuals summation : sum((model(x,xdata)-ydata).^2) . It is a number.
r: residuals only : model(x,xdata)-ydata .It is a vector variable
ef: exitflag which tells us the reason why the code stops 
    important numbers: 1: model is converged , everything is good.
                       0: number of iterations exceeds the limit
                       All other numbers are not good...
oput: output information of the least squared method results . It is a struct variable.
L: Lagrange multipliers at the boundaries... . It is also a struct variable.
jac: Jacobian : derivative values of the data points

Inorder to calculate chi2 and errorbars we have to code several lines here:
>> hes=jac'*jac; % This is the matrix where second derivatives are calculated.
>> variance_covariance=inv(hes); %inverse of covariance matrix
>> Cii=diag(variance_covariance); %diagonal elements of inv covariance matrix will be the
                                  %up and down acceptable values
>> noffreedom=length(X)-length(nuve); % number of freedom : N-m 
>> error_bar=full((1/1).*sqrt(((Cii)).*rnorm/(length(X)))*tinv(0.84135,noffreedom));%full converts sparse matrix attribute to regular double.
>> [
>> Afit1=gaussfit(nfit1,X);
we can also now check if rnorm is calculated correctly as it was explained above.
>> sum((Y-Afit1).^2) -> This gives the exact rnorm number. 

fprintf('\n Chi Square:%g ',chi_square);
fprintf('\n Amplitude:%g-+%g',nfit(1),nonzeros(error_bar(1)));
fprintf('\n Centre:%g-+%g',nfit(2),nonzeros(error_bar(2)));
fprintf('\n Width:%g-+%g',nfit(3),nonzeros(error_bar(3)));[
fprintf('\n Background:%g-+%g',nfit(4),nonzeros(error_bar(4)));                 

So what did happen ? Why did we have so small errorbars and chi2 numbers ?


One more thing: we did not use errorbar weighting in the fitting ?
How can we include it ? Well,  this infact is nt possible for many users but in reality we can still fit the data if we deal with a binamial data distribution where statistically the errorbar is same for all data points as sqrt of them.
Then we weight the data with the errorbar as ;
>> n1=[10 46 1 1];
>> dY=sqrt(Y); % only when the all standard deviations in the model or more preciesly in the data
% distribution is same... as it is in Gaussian, Lorenzian or in general Poisson Distribution datas.
% knowing this, we comfirm that , neutron scattering is a good example for that kind of distribution % where all errorbars are simply dA=sqrt(A).
% So having said this, we cannot do this trick for the errorbar set as it was before, when we % % %%%% created it as dY=rand(100,1).*2; % 
>> Yn=Y./dY;
>> [nfit2,rnorm2,r2,ef2,oput2,L2,jac2]=lsqcurvefit(@gaussfit,n1,X,Yn,lowerb,upperb,options);
>> Afit2=gaussfit(nfit2,X);
>> figure(11); hold on;
>> errorbar(X,Y,dY,'ob');
>> plot(X,Afit1,'-r'); % This the line where we do not have any errorbar
>> plot(X,Afit2.^2,'-g'); % This is the line where we have fit Yn but plot the square of it...
This is the trick where we use the weighting of the data ....
But the down side is the fit values are just for the sqrt(data)...
We still can calculate the real (true) parameter values :
1) nfit2_true(1)=nfit2(1)^2; %amplutudes 
2) nfit2_true(3)=nfit2(3)/sqrt(2); the width is needed to be devided by sqrt(2)
3) nfit2_true(2)=nfit2(2); all other values of fit parameters are same
4) nfit2_true(4)=nfit2(4); all other values of fit parameters are same
Again these corrections are due to the fact that we use the Amplitude Gaussian... If we had
used Area Gaussian we might have need to devide or square with different numbers.

>> hes2=jac2'*jac2;
>> variance_covariance2=inv(hes2);
>> Cii2=diag(variance_covariance2);
>> noffreedom2=length(X)-length(n1);
>>  error_bar2=full((1/1).*sqrt(((Cii2)).*rnorm2/(length(X)))*tinv(0.84135,noffreedom2))

 
___________________ REAL DATA _________________________________________________________
[verim]=specokuma('lc4O8_d0025_a',25); % the difficult data for a Gaussian Model 
q=verim(:,1);A=verim(:,end);dA=sqrt(verim(:,end)); 
figure(3); hold on;
errorbar(q,A,dA,'ob');

lowerb= [] ;
upperb=[];
options=optimset('Display','iter','MaxFunEvals',10000,'MaxIter',5000,'TolFun',1*10^(-4),'Algorithm','levenberg-marquardt');
n6=[300 0.22 0.001 100];
[nfit6,rnorm6,r6,ef6,oput6,L6,jac6]=lsqcurvefit(@gaussfit,n6,q,dA,lowerb,upperb,options);
% becareful here, we fit the data for dA... infact it is A./dA whichis simply dA. Right !!!!
Afit2dA=gaussfit(nfit6,q); 
plot(q,Afit2dA.^2,'-r');
% here we create our fit for the square of the original fit.

%now we can compare with non-LM way of fitting,
options=optimset('Display','iter','MaxFunEvals',10000,'MaxIter',5000,'TolFun',1*10^(-4));  
[nfit6,rnorm6,r6,ef6,oput6,L6,jac6]=lsqcurvefit(@gaussfit,n6,q,dA,lowerb,upperb,options);
% this fit will not be good as it was before...
Afit2dA=gaussfit(nfit6,q); 
plot(q,Afit2dA.^2,'-r');

chi_square=(rnorm6)/(noffreedom2)

