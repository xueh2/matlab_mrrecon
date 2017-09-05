
function data_list = PerformGadgetronRecon_FindPatient_From_SHA1(dataDir, sha1, study_date)
% data_list = PerformGadgetronRecon_FindPatient_From_SHA1(dataDir, sha1, study_date)
% data_list = PerformGadgetronRecon_FindPatient_From_SHA1('\\hl-share\RawMRI\Lab-Kellman\RawData\BARTS', '62085760c1c4de5dcd9926d569b39eee770d3b81', '20170215')
% given a sha1, find all data related to it

studyDir = fullfile(dataDir, study_date);

[names, num] = findFILE(studyDir, '*.h5');

data_list = [];

for ii=1:num
    try
        
        if( ~isempty(strfind(names{ii}, 'ISMRMRD_Noise_dependency')))
            continue;
        end
        
        dset = ismrmrd.Dataset(names{ii});
        header = ismrmrd.xml.deserialize(dset.readxml());
        id = header.subjectInformation.patientID;
        dset.close();
        
        disp([names{ii} ' - ' id]);
        
        if( ~isempty(strfind(id, sha1)) )        
            finfo = dir(names{ii});
            [path, name, ext] = fileparts(names{ii});
            
            data_list = [data_list; {name, finfo.bytes/1024/1024, study_date, names{ii} }];
        end        
    catch
        finfo = dir(names{ii});
        disp([names{ii} ' file size ' num2str(finfo.bytes/1024) 'k -- do not have patient ID ... '])
    end
end
