function [images, digits, angles] = digitTrain4DArrayData
% digitTrain4DArrayData   Load the digit training set as 4-D array data
%
%   images      - Input data as an H-by-W-by-C-by-N array, where H is the
%                 height and W is the width of the images, C is the number
%                 of channels, and N is the number of images.
%   digits      - Categorical vector containing the labels for each
%                 observation.
%   angles      - Numeric vector containing the angle of rotation in
%                 degrees for each image.

% Copyright 2016 The MathWorks, Inc.

[images, digits, angles] = digitTableToArray(digitTrainTable);

end