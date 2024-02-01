% In this trial, we are going to use GoogLeNet CNN for the library of the images.
% We are not going to train our CNN. 
% The steps are : 1) Load Camera Image & CNN 2) Resize the image for the Layer of the CNN
% Process the prediction. 4) Make [1:3] in a loop. 
% Written by Mmt on Feb. 2022 at ITU. 
 
kamera = webcam; % this is needed to Add-on from : USB WEBCAMS
net = googlenet; % This is the Deep Learning GoogLeNet Add-on 

boy=net.Layers(1).InputSize(1:2); % inputSize must be boy= [224   224];
figure;
resim=snapshot(kamera);
image(resim); % we can also use imshow(resim) ;
resimc=imresize(resim,boy); % now the original resim is cut!
[bul,basari]=classify(net,resimc);
title({char(bul),num2str(max(basari),2)});


kamera = webcam; % this is needed to Add-on from : USB WEBCAMS
net = googlenet; % This is the Deep Learning GoogLeNet Add-on 

ilk=-1e+2;
while ilk<10
ilk=ilk+1;
boy=net.Layers(1).InputSize(1:2); % inputSize must be boy= [224   224];
figure(10);
resim=snapshot(kamera);
%image(resim); % we can also use imshow(resim) ;
resimc=imresize(resim,boy); % now the original resim is cut!
image(resimc); % We did not cut by we degreased the resolution. 
drawnow; % this functions like hold on;
[bul,basari]=classify(net,resimc);
title({char(bul),num2str(max(basari),2)});
pause(1);
end;

