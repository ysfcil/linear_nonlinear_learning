function [images, digits, angles] = digitTableToArray( aDigitTable )
% digitTableToArray   Convert a table containing a digit dataset into a 4-D
% array
%
% Inputs
%   aDigitTable - A table containing training or test digit dataset.
%
% Outputs
%   images      - Input data as an H-by-W-by-C-by-N array, where H is the
%                 height and W is the width of the images, C is the number of
%                 channels, and N is the number of images.
%   digits      - Categorical vector containing the label of the digit for each
%                 observation.
%   angles      - Numeric vector containing the angle of rotation in
%                 degrees for each image.

% Copyright 2016 The MathWorks, Inc.

imagePaths = iGetFullDigitPaths( aDigitTable );

digitImds = imageDatastore( imagePaths, 'LabelSource', 'foldernames' );

[images, digits] = imds2array(digitImds);

angles = aDigitTable.angle;

end

function imagePaths = iGetFullDigitPaths( aDigitTable )
% Base digit path
digitPath = fullfile( matlabroot, 'toolbox', 'nnet', 'nndemos', 'nndatasets', 'DigitDataset' );

% Add digit folder to path name
imagePaths = strcat(num2str(aDigitTable.digit), filesep, aDigitTable.image);

% Add full path
imagePaths = cellfun(@(s)strcat(digitPath,filesep,s),imagePaths,'UniformOutput',false);
end
