%-------- INTORDUCTION TO  nlfit0  and ITS USAGE =================
[verim]=specokuma('lc4O8_d0025_a',29); 
q=verim(:,1);A=verim(:,end);dA=sqrt(verim(:,end));
figure(9); hold on;
errorbar(q,A,dA,'ob');

% the starting parameters : nuve % %%% as usual !!!!
n2=[3000 0.22 0.01 100]; 
mask=[1 1 1 1]; 

% the mask is new. We did not have this or similar neither for nlitfit nor lsqcurvefit .
% with this option, we fix any of the fitting parameter to the nuve values or let it
% recalculate and fit in every itterations....
 
[nfit,chisquare,errors,fitresult] = nlfit00_class(q,A,'gaussfit',n2,mask,dA);
% as you see, this command excepts specifically an errorbar matrix. IN other words, it is 
% completly a scientif tool, not like a biology gibrish  things !!!!! 

Afit=fitresult.yfit;   % it contains in a struct form everything we need, even the yfit values 
plot(q,Afit,'-r');     % as well.
                       % That means you do not need to recall the model with the fit-parameters
                       % one more time to calculate the yfit. It is done already for us to use it.

nuvefit
1416.6      0.22691     0.022811         1664
Chi-Square 8.047289

%Instancewise, if we fix the Amplitude in the gaussfit by putting "1" in the mask matrix on the
%exact position for amplitude;

mask=[0 1 1 1];
[nuvefit,chisquare,errors,fitresult] = nlfit00(q,A,'gaussfit',n2,mask,dA);
Afit=fitresult.yfit;
plot(q,Afit,'-g');

nuvefit
 3000 0.22331  0.007273  1861.2
Chi-Square 98.126320

% now let's try our new code with a difficult data to fit. We kknow it , it is # 25.
[verim1]=specokuma('lc4O8_d0025_a',25); 
q1=verim1(:,1);A1=verim1(:,end);dA1=sqrt(verim1(:,end));
figure(11); hold on;
errorbar(q1,A1,dA1,'sb');
nuve=[3000 0.22 0.01 100]; % let's use the same starting points... 
mask=[1 1 1 1];
[nfit1,chisquare1,errors1,fitresult1] = nlfit00(q1,A1,'gaussfit',nuve,mask,dA1); 

Afit1=fitresult1.yfit;
plot(q1,Afit1,'-r');

% so it is better and more robust then the command s of matlab that we had tried before....




