function y=gauss2D(n,M)
% 3D Gaussian Function Creator & Fitter
% written by MMT for FIZ508E Class Note in 2017


A=n(1);
w=n(2);
xc=n(3);
yc=n(4);

y=(A/w*sqrt(pi/2))*(exp(-2*((M(:,1)-xc).^2+(M(:,2)-yc).^2)./w^2));
