function noise_mear_id = findNoiseDependencyMeasurementID_SavedIsmrmrd(name)
% find the noise dependency ID of a saved ismrmrd scan
% noise_mear_id = findNoiseDependencyMeasurementID_SavedIsmrmrd(name)

try
    dset = ismrmrd.Dataset(name);
    header = ismrmrd.xml.deserialize(dset.readxml());

    noiseID = [];
    numDep = numel(header.measurementInformation.measurementDependency);
    for k=1:numDep
        if( strfind(header.measurementInformation.measurementDependency(k).dependencyType, 'Noise') == 1 )
            noiseID = header.measurementInformation.measurementDependency(k).measurementID;
            break;
        end
    end

    noise_mear_id = [];

    if(~isempty(noiseID))
        ind = strfind(header.measurementInformation.measurementID, '_');

        meas_id = header.measurementInformation.measurementID(1:ind(end)-1);
        if(numel(noiseID)<numel(meas_id))
            noise_mear_id = [meas_id '_' noiseID];
        else
            noise_mear_id = noiseID;
        end
    end
    
    dset.close();
catch
    noise_mear_id = [];
end