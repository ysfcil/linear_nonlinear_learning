kam=webcam; % this is done after adding MATLAB SUPPORT PACKAGE FOR USB WEBCAMSdeo
net=alexnet; %this is done after adding Deep Learning Toolbox Model For AlexNet Network

while 1 
resim=snapshot(kam);
image(resim);
resim1=imresize(resim,[227 227]);
cisim=classify(net,resim1);
fprintf('\n cisim=%s',cisim);
title(char(cisim));
drawnow; % this updates the figure continously and immediatly...
end;

% however in order to run a video -continous stream- from the camera
preview(kam);
