function timeUsed = run_gt_recon_case(dataName, configName, UTCases, deleteh5, startRemoteGT, res_suffix, h5Only, remoteXml, compressionBit, paraXml, UTDir)
% run the gt recon for a case
% timeUsed = run_gt_recon_case(dataName, configName, UTCases, deleteh5, startRemoteGT, res_suffix, h5Only, remoteXml, compressionBit, paraXml, UTDir)

if nargin<6
    res_suffix = [];
end

if nargin<7
    h5Only = 0;
end

if nargin<8
    remoteXml = 0;
end

if nargin<9
    compressionBit = 0;
end

if nargin<10
    paraXml = [];
end

if nargin<11
    UTDir = getenv('GTPLUS_UT_DIR')
end

OutputFormat = getenv('OutputFormat');

num = size(UTCases, 1);

for ii=1:num    
    ind = strfind(dataName, UTCases{ii, 2});
    if ( ~isempty(ind) & strcmp(configName, UTCases{ii, 4}) )
        testInfo = UTCases(ii, :);
        folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
        if(~isFileExist(folderDir))
            folderDir = fullfile(UTDir, UTCases{ii, 1})
        end
    
        if(~isFileExist(folderDir))
            folderDir = fullfile( UTCases{ii, 1}, UTCases{ii,2 });
        end
        
        dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
        
        if(~isempty(res_suffix))
            if(~h5Only)
                h5Name = fullfile(folderDir, [UTCases{ii, 2} res_suffix '.h5']);
            else
                h5Name = fullfile(folderDir, [UTCases{ii, 2} '.h5']);
            end
            resDir = fullfile(folderDir, [UTCases{ii, 5} res_suffix]);
        else
            h5Name = fullfile(folderDir, [UTCases{ii, 2} '.h5']);
            resDir = fullfile(folderDir, [UTCases{ii, 5}]);
        end
        
        xslFile = UTCases{ii, 7};

        isVD = strcmp(UTCases{ii, 3}, 'VD');
        isVD11 = strcmp(UTCases{ii, 3}, 'VD11');
        isAdjScan = strcmp(UTCases{ii, 3}, 'Adj');
        isNX = strcmp(UTCases{ii, 3}, 'NX');
        isNX20 = strcmp(UTCases{ii, 3}, 'NX20');
        
        delete(fullfile(folderDir, 'res*.h5'));
        delete(fullfile(folderDir, 'out*.h5'));
        delete(fullfile(folderDir, '*.xml'));

        if(strcmp(OutputFormat, 'hdr')==1)
            delete(fullfile(resDir, '*.nii'));
            delete(fullfile(resDir, 'gadgetron_*.hdr'));
            delete(fullfile(resDir, 'gadgetron_*.img'));
            delete(fullfile(resDir, 'Generic*.hdr'));
            delete(fullfile(resDir, 'Generic*.img'));
            delete(fullfile(resDir, 'GTPrep*.hdr'));
            delete(fullfile(resDir, 'GTPrep*.img'));
            delete(fullfile(resDir, 'GT*.hdr'));
            delete(fullfile(resDir, 'GT*.img'));
            delete(fullfile(resDir, 'Cmr*.hdr'));
            delete(fullfile(resDir, 'Cmr*.img'));
            delete(fullfile(resDir, '*.attrib'));
        end
        
        if(deleteh5)
            delete(fullfile(resDir, '*.xml'));
            delete(fullfile(resDir, '*.h5'));
        end
        
        timeUsed = run_gt_recon(folderDir, dataName, h5Name, deleteh5, isVD, isVD11, isNX, isNX20, isAdjScan, configName, resDir, xslFile, startRemoteGT, h5Only, remoteXml, compressionBit, paraXml);
    end    
end
