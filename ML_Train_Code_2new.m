% Written by Mmt Dec 2023 at ITU . It uses the preinstalled digitTrain4Darray set ! 
% It creates a linear DAG CNN for 15 layers. 

[resimler,labellar] = digitTrain4DArrayData; % Here imTrain is 28x28x1x5000 data.
                                         % All images are 28x28 resolution. There is 1 layer.
                                         % There are 5000 imaages !!! 
                                         % imValidates are the labels of each images . 

%imshow(resimler(:,:,:,12)); % creates the image for 9 
%labellar(12) % gives the categorical 9 ! 
                                         
idx = randperm(size(resimler,4),1000);  
                        %means rand(5000,1000);                                                  
                                         % size(XTrain,4) is 4000 , so get 1000 number randomly
                                         % from 4000 number, max. will 4000. 
resimkontrol = resimler(:,:,:,idx); % is now 28x28x1x1000 image file.
labelkontrol = labellar (idx);               % 1000 labels will be used to test the net 
%imshow(resimkontrol(:,:,1));              % draws the figure !
resimtrain=resimler;
resimtrain(:,:,:,idx)=[];
%resimler(:,:,:,idx) = [];                  %now these XTrain (1000) image is cleared! 28x28x1x1000
labeltrain=labellar;
labeltrain(idx)=[];
%labellar (idx) = [];                        %so labellar: 4000x1
                                         %now XTrain is 28x28x1x4000 & Ytrain is 4000 long!!

% We see that XValidation and YValidation are complaiant!

imageAugmenter = imageDataAugmenter('RandRotation',[-20,20], 'RandXTranslation',[-3 3], ...
    'RandYTranslation',[-3 3]);
%imageSize = [28 28 1];
imageSize=size(resimtrain(:,:,:,1));  % here size is 28 x 28 is alco working !!
augimds = augmentedImageDatastore(imageSize,resimtrain,labeltrain,'DataAugmentation',imageAugmenter);


layers = [imageInputLayer(imageSize)
          convolution2dLayer(3,8,'Padding','same')
          batchNormalizationLayer
          reluLayer   
          maxPooling2dLayer(2,'Stride',2)
          convolution2dLayer(3,16,'Padding','same')
          batchNormalizationLayer
          reluLayer   
          maxPooling2dLayer(2,'Stride',2)
          convolution2dLayer(3,32,'Padding','same')
          batchNormalizationLayer
          reluLayer   
          fullyConnectedLayer(10)
          softmaxLayer
          classificationLayer];

opts = trainingOptions("sgdm","ExecutionEnvironment","gpu","MaxEpochs",15, "Shuffle","every-epoch", ...
                       "Plots","training-progress","Verbose",false, ...
                        "ValidationData",{resimkontrol,labelkontrol});
  
%opts = trainingOptions("sgdm","ExecutionEnvironment","auto","MaxEpochs",15, "Shuffle",...
%                        "every-epoch", ...
%                       "Plots","training-progress","Verbose",false, ...
%                        "ValidationData",{resimkontrol,labelkontrol});

%opts = trainingOptions("sgdm","ExecutionEnvironment","cpu","MaxEpochs",15, "Shuffle",...
%                        "every-epoch", ...
%                       "Plots","training-progress","Verbose",false, ...
%                        "ValidationData",{resimkontrol,labelkontrol});

%opts = trainingOptions("sgdm","ExecutionEnvironment","parallel","MaxEpochs",15, "Shuffle",...
%                        "every-epoch", ...
%                         "Plots","training-progress","Verbose",false, ...
%                        "ValidationData",{resimkontrol,labelkontrol});

% delete(gcp('nocreate')); % To TURN OFF Multiple CPU option
% now we did not change anyting though ;;; 
%opts.ValidationData{1}(:,:,:,1)-XValidation(:,:,:,1) 

% gives only 0 's ... So they are same matrices. 
  
                       
net = trainNetwork(augimds,layers,opts);      

%parpool;
%tic;
%parfor i=1:gpuDeviceCount("available")
%net = trainNetwork(augimds,layers,opts);
%end
%toc
augimdsValidation = augmentedImageDatastore(imageSize,resimkontrol);
classify(net,resimkontrol(:,:,3)) % now we can use the net to define what have in.
[Predictions, Probabilities]= classify(net,augimdsValidation);
accuracy = mean(Predictions == labelkontrol); 
%[YPred,probs] = classify(net,augimds);  
%accuracy = mean(YPred == labelkontrol); % gives out 0.9675 
