function openh5(filename)
% Wrapper function to open HDF5 files in SimpleViewer when
% drag-and-dropped into MATLAB
%
% Kelvin Chow (kelvin.chow@virginia.edu)
% Department of Cardiovascular Medicine
% University of Virginia Health System
% Revision: 1.0.1  Date: 16 January 2015
%
% Copyright (c) 2015 Kelvin Chow

	figure;
	SimpleViewer(filename);
end