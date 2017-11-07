
% idir = 'Z:\RawData\RoyalFree\'
% startdate = 20161015;
% enddate = 20170724;
idir = '\\hl-share\RawMRI\Lab-Kellman\RawData\LEEDS\';
startdate = 20170601;
enddate = 20171026;

clear out
        
k = 1; m = 1;
d1 = dir([idir,filesep,'20*']);
for i = 1:length(d1)
    dirname = d1(i).name;
    date = str2num(dirname);
    if date >= startdate & date <= enddate
        
        dirname = [idir,filesep,d1(i).name,filesep];
        
        db_psir_list = dir([dirname,'*Perf*h5']);
        for j = 1:length(db_psir_list)
            filename = db_psir_list(j).name;
%             disp(filename)
            s = dir([dirname,filename]);
            if s.bytes > 0; % 1e9 % sax stack
                [patientID, protocol] = read_ismrmrd_protocol([dirname,filename]);
                if ~isempty(protocol)
%                     disp([d1(i).name, '    ',tmp.seqname])
                    
                    ind = strfind(filename, '_');
                    str1 = filename(ind(5)+1:ind(8)-1);
                    ind = strfind(filename, '-');
                    studydir{k} = dirname;
                    disp([studydir{k},' ',filename,'  ', num2str(protocol.lRepetitions+1)])
                    out(k).studydir = studydir{k};
                    out(k).filename = filename;
                    out(k).measurements = protocol.lRepetitions+1;
                    k = k + 1;
                end
            end
        end
    end
end

return
                % find perfusion
%                 p = dir ([dirname,filesep,'Perf*',str1,'*']); 
%                 if ~isempty(p);
%                     perf_filename = p(end).name;
%                     ind = strfind(perf_filename, '-');
%                     perf_time(k) = hhmmss2time(perf_filename(ind(1)+1:end-3));
%                     if perf_time(k) > db_time(k) % must have been a retro recon at end of the study
%                         perf_filename = p(end-1).name;
%                         ind = strfind(perf_filename, '-');
%                         perf_time(k) = hhmmss2time(perf_filename(ind(1)+1:end-3));                        
%                     end
%                 else
%                     perf_time(k)=0;
%                 end
%                 bb_lge = dir ([dirname,filesep,'LGE*',str1,'*']); 
%                 if ~isempty(bb_lge);
%                     clear stack
%                     % find sax stack
%                     for n = 1:length(bb_lge)
%                         stack(n) = bb_lge(n).bytes > 1e9;
%                     end
%                     m = find(stack==1);
%                     if ~isempty(m)
%                         lge_filename = bb_lge(m(end)).name;
%                         ind = strfind(lge_filename, '-');
%                         lge_time(k) = hhmmss2time(lge_filename(ind(1)+1:end-3));
%                     else
%                         lge_time(k)=0;
%                     end
%                 else
%                     lge_time(k)=0;
%                 end                
%                                 
%                 k = k + 1;
%             end
%             m = m + 1;
%         end
%     end
% end


% clear te td1 td2 delta t1myo t1blood
% for k = 1:length(tmp)
%     te(k) = tmp{k}.TE;
%     td1(k) = tmp{k}.TD1;
%     td2(k) = tmp{k}.TD2;
%     delta(k) = tmp{k}.delta;
%     t1myo(k)= tmp{k}.T1_myo;
%     t1blood(k)= tmp{k}.T1_blood;
% end
% 
% disp(['N  = ',num2str(length(te))])
% disp(['TE  = ',num2str(mean(te)),' +/- ',num2str(std(te))])
% disp(['TD1 = ',num2str(mean(td1)),' +/- ',num2str(std(td1))])
% disp(['TD2 = ',num2str(mean(td2)),' +/- ',num2str(std(td2))])
% disp(['T1myo = ',num2str(mean(t1myo)),' +/- ',num2str(std(t1myo))])
% disp(['T1blood = ',num2str(mean(t1blood)),' +/- ',num2str(std(t1blood))])
% 

% 
% clear lge db
% for k = 1:length(db_time)
%     if perf_time(k) > 0 & lge_time(k) >0
%         lge(k) = lge_time(k) - perf_time(k);
%         db(k) = db_time(k) - perf_time(k);
%     else
%         lge(k) = 0; db(k)=0;
%     end
% end
% 
% db = db/60; lge = lge/60;
% 
% figure; subplot(2,1,1);hist(db(db>0),32); title ('DB PSIR')
% subplot(2,1,2);hist(lge(db>0),32); title (' Bright blood PSIR')
% 
% 
% 
% t_db = db(db>0);
% t_bb = lge(db>0);
% 
% bb_1st = find(t_bb < t_db);
% db_1st = find(t_bb >= t_db);
% 
% 
% m1 = mean(t_bb(db_1st));
% m2 = mean(t_db(db_1st));
% m3 = mean(t_bb(bb_1st));
% m4 = mean(t_db(bb_1st));
% 
% s1 = std(t_bb(db_1st));
% s2 = std(t_db(db_1st));
% s3 = std(t_bb(bb_1st));
% s4 = std(t_db(bb_1st));
% 
% disp(['DB first:'])
% disp(['    BB: ',num2str(m1),' +/- ',num2str(s1)])
% disp(['    DB: ',num2str(m2),' +/- ',num2str(s2)])
% disp(['BB first:'])
% disp(['    BB: ',num2str(m3),' +/- ',num2str(s3)])
% disp(['    DB: ',num2str(m4),' +/- ',num2str(s4)])





