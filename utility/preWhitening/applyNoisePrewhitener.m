function outputdata = applyNoisePrewhitener(inputdata, noisePrewhiteningMatrix);
% function outputdata = applyNoisePrewhitener(inputdata, noisePrewhiteningMatrix);
%
% function to apply noise pre-whitener matrix
%
% inputs:
%     inputdata assumes that Channels are 1st dimension
%     noisePrewhiteningMatrix is nChannels x nChannels matrix
% outputs:
%     outputdata is noise pre-whitened arranged with Channels as 1st dimension

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

s=size(inputdata);
inputdata  = reshape(inputdata, [s(1) prod(s(2:end))]);
outputdata = noisePrewhiteningMatrix * inputdata;
outputdata = reshape(outputdata, s);

return





