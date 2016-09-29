function [data, header, readouts] = ismrmrd_read_data(datafile)
% open and read in ismrmrd h5 data file

data = ismrmrd.Dataset(datafile);
header = ismrmrd.xml.deserialize(data.readxml());

if(nargout>2)
    readouts = data.readAcquisition();
end