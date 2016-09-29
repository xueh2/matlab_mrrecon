
function PerformGadgetronRecon_SendDicom(resDir, date_start, date_end, dicomServer, dicomPort)
% PerformGadgetronRecon_SendDicom(resDir, date_start, date_end, dicomServer)
% PerformGadgetronRecon_SendDicom('I:\ReconResults\BARTS', '2016-05-12', '2016-06-14', 'barbados')

if(nargin<2)
    date_start = '2016-01-01';
end

if(nargin<3)
    date_end = '2017-01-01';
end

if(nargin<4)
    dicomServer = 'barbados';
end

if(nargin<5)
    dicomPort = 11112;
end

% ------------------------------------------------------------

startN = datenum(date_start);
endN = datenum(date_end);

[subdirs, num] = FindSubDirs(resDir);

for ii=1:num
    
    tt = datenum(subdirs{ii}, 'yyyymmdd');  
    if (tt<=endN && tt>=startN)
        [subdirs_dicom, num_dicom] = FindSubDirs(fullfile(resDir, subdirs{ii}));
        
        for kk=1:num_dicom
            
            if ( ~isempty(strfind(subdirs_dicom{kk}, 'dicom')) )
            
                dcmFolder = fullfile(resDir, subdirs{ii}, subdirs_dicom{kk});
                disp(['Sending ' dcmFolder ' to ' dicomServer]);
                
                command = ['D:\gtuser\gt_windows_setup\dcmtk-3.6.0\install_vc14\bin\storescu ' dicomServer '.nhlbi.nih.gov ' num2str(dicomPort) ' ' dcmFolder ' --scan-directories -aec DCM4CHEE --user admin --password admin'];
                tic; dos(command, '-echo'); toc
            end
        end
    end    
end
