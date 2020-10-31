function scanner_list = print_saved_Ismrmrd_scanner_id(flag_print);

% print_saved_Ismrmrd_scanner_id

if(nargin<1)
    flag_print = 1;
end

sites = {'BARTS'; 'ROYALFREE'; 'KCL'; 'LEEDS'; 'LUND'; 'KAROLINSKA'; 'CHENIESMEWS'; 'CAMPINAS'; 'GEISINGER'; 'CATHLAB'; 'CNMC'};
scanner_id = { [41837, 42110, 66016]; 42363; 42170; 66097; 42034; [41672, 46184]; [141303, 166032]; 53531; 42311; 41624; 41548};

scanner_list =  table(sites, scanner_id);

if(flag_print)
    for ii=1:size(scanner_list)
        disp([scanner_list.sites{ii} ' --------- ' num2str(scanner_list.scanner_id{ii})]);
    end
end