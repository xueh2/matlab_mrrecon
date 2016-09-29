function timeInSeconds = ConvertDicomAcquisitionTime2Seconds(cAcquisitionTime)
% bool convertAcquisitionTime2Second(const char* cAcquisitionTime, 
%                                      int& iHours, int& iMinutes, 
%                                      int& iSeconds, double& dSecRes,
%                                      double& timeInSeconds)
% {
%     try
%     {
%         double dAcquisitionTime;
% 
%         char *cStopString;       
%         dAcquisitionTime = strtod(cAcquisitionTime,&cStopString);
%         
%         if (!strncmp(cAcquisitionTime,cStopString, sizeof (cAcquisitionTime)))
%         {
%             FTK_ERROR_MESSAGE_STREAM("ERROR: conversion of Acquisition Time failed in FMRIFunctor::convertRelTime_ms():"<<cAcquisitionTime);
%             return false;
%         }
%         
%         // DICOM 24h-time format: hhmmss.ssssss
%         iHours   = int(  dAcquisitionTime / 10000. );
%         iMinutes = int( (dAcquisitionTime - iHours*10000.) / 100. );
%         iSeconds = int( (dAcquisitionTime - iHours*10000. - iMinutes*100.) );
%         dSecRes  = dAcquisitionTime - iHours*10000. - iMinutes*100. - iSeconds;
% 
%         timeInSeconds = iHours*3600 + iMinutes*60 + iSeconds + dSecRes;
%     }
%     catch(...)
%     {
%         FTK_ERROR_MESSAGE_STREAM("Exceptions happened in convertAcquisitionTime2Second(...) ... ");
%         return false;
%     }
% 
%     return true;
% }

ind = find(cAcquisitionTime=='.');

hhmmss = cAcquisitionTime(1:ind-1);
secStr = cAcquisitionTime(ind+1:end);

hour = str2double(hhmmss(1:2));
minute = str2double(hhmmss(3:4));
second = str2double(hhmmss(5:6));
subSec = str2double(['0.' secStr]);

timeInSeconds = hour*3600 + minute*60 + second + subSec;
