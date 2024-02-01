kam=webcam;% this is done after adding MATLAB SUPPORT PACKAGE FOR USB WEBCAMS

videoyaz=VideoWriter('videodeneme1.avi');

open(videoyaz); % opens the avi file

tic;
for ilk=1:600
 resim=snapshot(kam);
 writeVideo(videoyaz,resim);
 image(resim);
 drawnow;
 fprintf('\n ilk=%g',ilk);
 pause(0.01);
end;
toc;  % this means 60 ilk makes 6 s video capturing.
close(videoyaz);

clear('kam');
 
 %but his video when it is played in player, it is only 2 second long.
