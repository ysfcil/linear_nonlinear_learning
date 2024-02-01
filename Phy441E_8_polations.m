% SHARED CLASS NOTES -- INTERPOLATION & EXTRAPOLATION -- 
% These two techniques; interpolation and exterpolations are very useful if they are used % % % % %  logically. This ststement is especially % % valid for the later one.  
% we will the effect of interpolation first on the results of a Gaussian Fit. Then, we will
% what would we obtain from exterpolation. 
% Now we need a less frequecied data (meaning a data with less data point) for this purpose : 

[verim]=specokuma('lc4O8_d0025_a',1) ; 
q=verim(:,1);A=verim(:,end);dA=sqrt(verim(:,end));

% the number of that data point here is:  21 !!!
% First we need to fit that the data without any further manupulation.  

nuve=[10000 1 0.1 100];
mask=[1 1 1 1];
[nfit,chisquare,errors,fitresult] = nlfit00(q,A,'gaussfit',nuve,mask,dA);

figure(1); hold on; 
errorbar(q,A,dA,'ko');
plot(q,Afit,'b-');

% here notice here two very important thing : 1) the curve is not smooth especially at the peak  
% peak position 2) Chi-square is 45 !!! 

% now let's interpolate the data before we fit it. Let's increase the frequency od sampling by 10

adim_size=(q(end)-q(1))/length(q);
qinterp=q(1):adim_size/10:q(end);
% or simply -> 
qinterp2=q(1):(mean(diff(q)))/10:q(end);
% here we created the x-axis sample positions not the y-axis values yet. 

Ainterp=interp1(q,A,qinterp2); %this uses default thecnique of linear connections
%or using : cubic or pchip technique to obtaine SMOOTH curves
Ainterp1=interp1(q,A,qinterp,'pchip');
Ainterp1=interp1(q,A,qinterp,'cubic');

dAinterp=sqrt(Ainterp); % now we have the interpolated errorbars of data. 
errorbar(qinterp2,Ainterp,dAinterp,'r.'); % the new data shown by red dots. 

% Now let's fit the data one more time bt this time using the interpolated data. 
[nfit2,chis2,errors2,fitresult2] = nlfit00(qinterp2,Ainterp,'gaussfit',nuve,mask,dAinterp);
Afit2=fitresult2.yfit;
plot(qinterp2,Afit2,'g-',linewidth',2);

%Now notice two things again : 1) The green line is much smoother than the blue one 
 % and 2) our new chi-square decreased to 31. 
   
% Again we could try something else : Like, we might interpolate the fitted data (original
% fit line by 10 times more and check the result... OK, Let's do it next ;
Afit1interp=interp1(q,Afit,qinterp2); % this increases the number of Afit (original) one by 10
plot(qinterp2,Afit1interp,'c-');

% Now we get again not a smooth line. The shape looks exactly like the original fit line with
% sharp peak lines. 
% So, what we learn that interpolation creates a mirage of better fit line and a better 
% chi-square numbers. If the person - Pepople do not know the details , you can simply 
% make them believe that what you are doing nice and stringly correct. 
% in fact the strictly correct fit is the first one, with sharp lines and higher chisquares.
 
% We might try spline method on interpolation .
 
Afit1a=interp1(q,Afit,qinterp2,'spline');
% This might solve your sharp line peak position problem but the figure looks a bit more
% ridiculous, and suspicious of having less chis-square. 


%%%%% ----- %%%%%%%%%%%%%%%%%%%%%%%%% ---- %%%%%%%%%%%%%
%%%% % EXTRAPOLATION % %%%%%%%%%%%%%%


step_size=(mean(diff(q)))/2;
qext=[qinterp,qinterp(end)+step_size,qinterp(end)+(step_size*2),qinterp(end)+(step_size*3)];
Aext=interp1(q,A,qext, 'linear', 'extrap');

figure(5); hold on;
plot(q,A,'ob');
plot(qext,Aext,'or');
