function [yfitmin]=fitAlan(n,X) 
%size(X);
%size(beta);

%pause;

bg=n(1);
Amp=abs(n(2));
xc=n(3);
wx=abs(n(4));
yc=n(5);
wy=abs(n(6));

%disp('simdi pause var!');pause;

%mx=reshape(X(:,1),250,150);%mx=reshape(X(:,1),150,100);
%my=reshape(X(:,2),250,150);%my=reshape(X(:,2),150,100);
%mz=reshape(y,150,100);
%[xd,yd]=meshgrid(xad,yad);
%xd=mx;yd=my;
zd=Amp*exp(-1*(((X(:,1)-xc).^2/wx^2)+((X(:,2)-yc).^2/wy^2)))+bg;
zd=zd(:); % we know that zd is ready, but for caution, we used one linearization line here.
%disp('puse var!');pause;
%yfitmin=sum(sum((abs(zd)-abs(A)).^2));
yfitmin=zd;
