
function res = PerformGadgetronRecon_parse_scan_info(scan_info)
% parse the scan info
% res = PerformGadgetronRecon_parse_scan_info(scan_info)

res = struct('num_pt', 0, 'LGE', 0, 'pt_LGE', 0, 'DBLGE', 0, 'pt_DBLGE', 0, 'Perf', 0, 'pt_Perf', 0, 'Binning', 0, 'pt_Binning', 0, 'RTCine', 0, 'pt_RTCine', 0, 'T2W', 0, 'pt_T2W', 0, 'T2S', 0, 'pt_T2S', 0, 'study_dates', [], 'week_num', []);

res.num_pt = size(scan_info, 1);

res.LGE = sum(scan_info.LGE);
res.pt_LGE = sum(scan_info.LGE>0);

res.DBLGE = sum(scan_info.DBLGE);
res.pt_DBLGE = sum(scan_info.DBLGE>0);

res.Perf = sum(scan_info.Perf);
res.pt_Perf = sum(scan_info.Perf>0);

res.Binning = sum(scan_info.Binning);
res.pt_Binning = sum(scan_info.Binning>0);

res.RTCine = sum(scan_info.RTCine);
res.pt_RTCine = sum(scan_info.RTCine>0);

res.T2W = sum(scan_info.T2W);
res.pt_T2W = sum(scan_info.T2W>0);

res.T2S = sum(scan_info.T2S);
res.pt_T2S = sum(scan_info.T2S>0);

start_date = [];
end_date = [];

mint = 0;
maxt = 0;

minw = 0;
maxw = 0;

for ii=1:res.num_pt 
    curr_date = scan_info.study_dates(ii);
    t = curr_date;
    
    w = scan_info.week_num(ii);
    
    if(ii==1)
        mint = t;
        maxt = t;
        start_date = curr_date;
        end_date = curr_date;
        
        minw = scan_info.week_num(ii);
        maxw = scan_info.week_num(ii);
    else
        if(t<mint)
            mint = t;
            start_date = curr_date;
        end
        if(t>maxt)
            maxt = t;
            end_date = curr_date;
        end
        
        if(w<minw)
            minw = w;
        end
        if(w>maxw)
            maxw = w;
        end
    end    
end

res.study_dates = [start_date end_date];
res.week_num = [minw maxw];

disp(['---------------------------------------------------']);
disp(['Total patients   : ' num2str(res.num_pt)]);
disp(['Time range       : ' num2str(start_date) ' to ' num2str(end_date)]);
disp(['LGE              : ' num2str(res.LGE) ' from ' num2str(res.pt_LGE) ' patients']);
disp(['DB LGE           : ' num2str(res.DBLGE) ' from ' num2str(res.pt_DBLGE) ' patients']);
disp(['Perfusion        : ' num2str(res.Perf) ' from ' num2str(res.pt_Perf) ' patients']);
disp(['Binning          : ' num2str(res.Binning) ' from ' num2str(res.pt_Binning) ' patients']);
disp(['RT Cine          : ' num2str(res.RTCine) ' from ' num2str(res.pt_RTCine) ' patients']);
disp(['T2W              : ' num2str(res.T2W) ' from ' num2str(res.pt_T2W) ' patients']);
disp(['T2S              : ' num2str(res.T2S) ' from ' num2str(res.pt_T2S) ' patients']);
disp(['---------------------------------------------------']);

