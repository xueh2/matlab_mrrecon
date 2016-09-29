
function [tUsed, ignored] = PerformGadgetronRecon_Fix_SavedIsmrmrd(dataDir, date_start, date_end)
% [tUsed, ignored] = PerformGadgetronRecon_Fix_SavedIsmrmrd(dataDir, date_start, date_end)
% [tUsed, ignored] = PerformGadgetronRecon_Fix_SavedIsmrmrd('I:\KAROLINSKA')
% [tUsed, ignored] = PerformGadgetronRecon_Fix_SavedIsmrmrd('I:\ROYALFREE.OLD')
% [tUsed, ignored] = PerformGadgetronRecon_Fix_SavedIsmrmrd('I:\BARTS')

date_suffix = datestr(date, 'yyyymmdd');
studyDate = datestr(date, 'yyyy-mm-dd');

if(nargin<2)
    date_start = '2016-01-01';
end

if(nargin<3)
    date_end = '2017-01-01';
end

% ------------------------------------------------------------

% find data

files = [];

configNames = [];
startN = datenum(date_start);
endN = datenum(date_end);

[names, num] = findFILE(dataDir, '*.h5');          
for n=1:num
           
    [pathstr, name, ext] = fileparts(names{n});
    
    % find scanner ID, patient ID, study ID, measurement ID, study date and time
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
       
    tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));
    
    if (tt<=endN && tt>=startN)
        disp(name);
        files = [files; {name}];
        configNames = [configNames; {configName}];
    end
end 

num = numel(files);
tUsed = [];
ignored = [];
for n=1:num

    name = files{n};   
    
    if( isempty(strfind(name, 'Noise')) == 0 )
        continue;
    end
    
    dataName = fullfile(dataDir, [name '.h5']);
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);             
    
    %% fix the scan

    try
        dset = ismrmrd.Dataset(dataName);
        hdr = ismrmrd.xml.deserialize(dset.readxml);

        D = dset.readAcquisition(1, 1);
    catch
        disp([name ' is read in wrongly ... '])
        % dset.close();
        continue;
    end
    timeInSeconds = double(D(1).head.acquisition_time_stamp) * 2.5/1e3;
    
    hours = floor(timeInSeconds/3600);
    mins =  floor((timeInSeconds - hours*3600) / 60);
    secs =  floor(timeInSeconds- hours*3600 - mins*60);
    
    fix_header = 0;
    fix_study_date = 0;
    fix_study_date_using_data = 0;
    fix_study_time = 0;
    fix_study_time_using_data = 0;
    studyDate_correct_format = 1;
    studyTime_correct_format = 1;
    
    if(~isfield(hdr, 'studyInformation'))
        fix_study_date = 1;
        fix_study_time = 1;
    else
        if(~isfield(hdr.studyInformation, 'studyDate'))
            fix_study_date = 1;
        else
            if(strcmp(hdr.studyInformation.studyDate, [study_year study_month study_day])==0 && strcmp(hdr.studyInformation.studyDate, [study_year '-' study_month '-' study_day])==0)
                fix_study_date = 1;
                fix_study_date_using_data = 1;
                
                ind_ = strfind(hdr.studyInformation.studyDate, '-');
                if(~isempty(ind_))
                    studyDate = [hdr.studyInformation.studyDate(1:4) hdr.studyInformation.studyDate(6:7) hdr.studyInformation.studyDate(9:end)];
                else
                    studyDate = hdr.studyInformation.studyDate;
                end
            else
                if(isempty(strfind(hdr.studyInformation.studyDate, '-')))
                    studyDate_correct_format = 0;
                else
                    studyDate_correct_format = 1;
                end
            end
        end
        
        if(~isfield(hdr.studyInformation, 'studyTime'))
            fix_study_time = 1;
        else
            if(strcmp(hdr.studyInformation.studyTime, [study_time(1:2) study_time(3:4) study_time(5:6)])==0 && strcmp(hdr.studyInformation.studyTime, [study_time(1:2) ':' study_time(3:4) ':' study_time(5:6)])==0)
                fix_study_time = 1;
                fix_study_time_using_data = 1;
               
                ind_ = strfind(hdr.studyInformation.studyTime, ':');
                if(~isempty(ind_))
                    studyTime = [hdr.studyInformation.studyTime(1:2) hdr.studyInformation.studyTime(4:5) hdr.studyInformation.studyTime(7:end)];
                else
                    studyTime = hdr.studyInformation.studyTime;
                end
            else
                if(isempty(strfind(hdr.studyInformation.studyTime, ':')))
                    studyTime_correct_format = 0;
                else
                    studyTime_correct_format = 1;
                end
            end
        end        
    end
      
    if(fix_study_date)
        if(~fix_study_date_using_data)
            hdr.studyInformation.studyDate = [study_dates(1:4) '-' study_dates(5:6) '-' study_dates(7:8)];
        end
    end
    
    if(fix_study_time)
        if(~fix_study_time_using_data) 
            hdr.studyInformation.studyTime = [num2str(hours, '%02i') ':' num2str(mins, '%02i') ':' num2str(secs, '%02i')];
        end
    end
    
    if(~isempty(strfind(name, 'Perf')))
        hdr.encoding(2).trajectoryDescription.identifier = 'AIF';
    end
    
    if(fix_study_date | fix_study_time)
       
        indA = strfind(name, '_');
        indB = strfind(name, '-');
        if(fix_study_date_using_data)            
            new_dataName = [name(1:indA(end)) studyDate '-'];
        else
            new_dataName = name(1:indB(end));
        end

        if(fix_study_time_using_data) 
            new_dataName = [new_dataName studyTime '.h5'];
        else
            new_dataName = [new_dataName num2str(hours, '%02i') num2str(mins, '%02i') num2str(secs, '%02i') '.h5'];
        end
        
        fprintf(2,'%d out of %d - Fixed - %s --> %s\n', n, num, name, new_dataName);
        %disp([num2str(n) ' out of ' num2str(num) '- Fixed - ' name ' --> ' new_dataName]);
        
        if(studyDate_correct_format==0)
            sd = hdr.studyInformation.studyDate;
            hdr.studyInformation.studyDate = [sd(1:4) '-' sd(5:6) '-' sd(7:8)];
        end
        
        if(studyTime_correct_format==0)
            st = hdr.studyInformation.studyTime;
            hdr.studyInformation.studyTime = [st(1:2) ':' st(3:4) ':' st(5:6)];
        end
        
        xmlstring = ismrmrd.xml.serialize(hdr);
        dset.writexml(xmlstring);
    else
        if(studyDate_correct_format==0 | studyTime_correct_format==0)
            
            if(studyDate_correct_format==0)
                sd = hdr.studyInformation.studyDate;
                hdr.studyInformation.studyDate = [sd(1:4) '-' sd(5:6) '-' sd(7:8)];
            end
            
            if(studyTime_correct_format==0)
                st = hdr.studyInformation.studyTime;
                hdr.studyInformation.studyTime = [st(1:2) ':' st(3:4) ':' st(5:6)];
            end
            
            xmlstring = ismrmrd.xml.serialize(hdr);
            dset.writexml(xmlstring);
        
            fprintf(2,'%d out of %d - Fixed - %s --> %s - %s\n', n, num, name, hdr.studyInformation.studyDate, hdr.studyInformation.studyTime);            
            % disp([num2str(n) ' out of ' num2str(num) '- Fixed - ' name ' --> ' hdr.studyInformation.studyDate ' - ' hdr.studyInformation.studyTime]);
        else
            disp([num2str(n) ' out of ' num2str(num) '- NOT Fixed - ' name]);
        end
    end
    
    dset.close();
    
    if(fix_study_date | fix_study_time)
        
        if(~isFileExist(fullfile(dataDir, new_dataName)))
            command = ['rename ' dataName ' ' new_dataName]
            dos(command, '-echo');
        else
            disp([num2str(n) ' out of ' num2str(num) '- already exists --> ' new_dataName]);
        end
    end
end
