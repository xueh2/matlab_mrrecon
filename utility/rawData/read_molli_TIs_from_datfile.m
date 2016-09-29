function TI = read_molli_TIs_from_datfile(datfilename);

% datfilename='D:\data\HUI\MOLLI RAW DATA\20101028_12h23m22s_54639.dat'

readout=mdh(datfilename); 
for i=1:length(readout)
    [header] = read_n4_header_ver4(datfilename, readout, i);
    aushFreePara(:,i)=header.aushFreePara;
end
evalInfoMask = readout(:,9);
MDH_LASTSCANINSLICE = find(bitand(evalInfoMask, 2^29)~=0);
TI = aushFreePara(1,MDH_LASTSCANINSLICE);
