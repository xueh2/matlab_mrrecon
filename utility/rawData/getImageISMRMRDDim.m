
function [slc, e2, con, phs, rep, set, ave, cha, run, imagenum, series] = getImageISMRMRDDim(name)
% get the ISMRMRD dimensions

indRUN = strfind(name, '_Run');
indSLC = strfind(name, '_SLC');
indE2 = strfind(name, '_E2');
indCON = strfind(name, '_CON');
indPHS = strfind(name, '_PHS');
indREP = strfind(name, '_REP');
indSET = strfind(name, '_SET');
indAVE = strfind(name, '_AVE'); indAVE = indAVE(end);
indCHA = strfind(name, '_CHA');
ind = strfind(name, '_');

if ( isempty(indRUN) )
    run = 0;
else   
    run = str2double(name(indRUN(end)+4:ind(2)-1));
end

slc = str2double(name(indSLC(end)+4:indCON(end)-1));
% e2  = str2double(name(indE2+3:indCON-1));
e2 = 0;
con = str2double(name(indCON(end)+4:indPHS(end)-1));
phs = str2double(name(indPHS(end)+4:indREP(end)-1));
rep = str2double(name(indREP(end)+4:indSET(end)-1));
set = str2double(name(indSET(end)+4:indAVE(end)-1));

N = numel(ind);

ave = str2double(name(indAVE(end)+4:ind(N-1)-1));
% cha = str2double(name(indCHA+4:indSeries-1));
cha = 0;

imagenum = str2double(name(ind(N-1)+1:ind(N)-1));
series = str2double(name( ind(N)+1:end));