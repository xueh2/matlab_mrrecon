function [readoutindicator]=mdh(datfilename);
% function readoutindicator=mdh(datfilename);
%
% function that plots the phase encode order vs readout # (time)
%
% usage:
%     readoutindicator=mdh(datfilename);
%
%     phaseencode     =  readoutindicator(:,1);
%     acquisition     =  readoutindicator(:,2);
%     slice           =  readoutindicator(:,3);
%     partition       =  readoutindicator(:,4);
%     echo            =  readoutindicator(:,5);
%     phase           =  readoutindicator(:,6);
%     repetition      =  readoutindicator(:,7);
%     dataset         =  readoutindicator(:,8);
%     ulEvalInfoMask  =  readoutindicator(:,9);
%     ulTimeStamp     =  readoutindicator(:,10);
%     ulPMUTimeStamp  =  readoutindicator(:,11);
%     readoutposition =  readoutindicator(:,12);
%     ushSamplesInScan=  readoutindicator(:,13);
%     TSE seq counter =  readoutindicator(:,16);
%     Ice Dimension a =  readoutindicator(:,17);
%     Ice Dimension b =  readoutindicator(:,18);
%     Ice Dimension c =  readoutindicator(:,19);
%     Ice Dimension d =  readoutindicator(:,20);
%     Ice Dimension e =  readoutindicator(:,21);


%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

[headers,protocol_header]=read_dat_headers(datfilename);

switch protocol_header.ulVersion
    case ' 0x3e9'
        version=15;
    case ' 0x7da'
        version=21;
    case ' 0xbee332'
        version=25;
    case ' 0x1421cf5' % vb11
        version=11;
    case ' 0x1452a3b' %vb13
        version=13;
    case ' 0x1483779' % vb15a
        version=13;
    case ' 0x1485e85' % vb15b
        version=13;
    case ' 0x14b44b6' % VB17A
        version=13;
end

[rawheader, readoutindicator, readoutindicator_full] = ...
    read_n4_indicators_ver3a(datfilename, version);

readoutindicator=readoutindicator_full;

return

figure; plot(readoutindicator(:,1),'.')
ylabel('phase encode')
xlabel('time')
title('phase encode order')

h_fig=figure;

h_axis(1)=subplot(1,3,1); 
plot(readoutindicator(:,1),'.')
ylabel('phase encode')
xlabel('time')

h_axis(2)=subplot(1,3,2);
plot(readoutindicator(:,4),'.')
ylabel('partition')
xlabel('time')

h_axis(3)=subplot(1,3,3); 
plot(readoutindicator(:,3),'.')
ylabel('slice')
xlabel('time')

linkaxes(h_axis,'x');
zoom xon
pan(h_fig,'xon');

[path,name]=fileparts(rawfilename);
name=strrep(name,'_',' ');
% title(name);


% add other plots
% set position/size
% add xaxis
% add figure name










