function trainTable = digitTrainTable
digitTrainCSVPath = fullfile( matlabroot, 'toolbox', 'nnet', 'nndemos', 'nndatasets', 'DigitDataset', 'digitTrain.csv' );
trainTable = readtable( digitTrainCSVPath, 'Delimiter', ',' );
end

