function [RampupTime, RampdownTime, FlattopTime, DelaySamplesTime, ADCDuration, DestSamples] = findRampSamplingTimes(header)
% function [RampupTime, RampdownTime, FlattopTime, DelaySamplesTime, ADCDuration, DestSamples] = findRampSamplingTimes(header);

RampupTime = 0;
RampdownTime = 0;
FlattopTime = 0;
DelaySamplesTime = 0;
ADCDuration = 0;
DestSamples = 0;

str = '<ParamLong."RampupTime">';
ind = strfind(header, str);
len = numel(str);
str2 = '<Visible> "true"';
ind2 = strfind(header(ind+len:ind+len+200), str2);

if ( isempty(ind2) )
    return;
end

startInd = ind+len+ind2(1)+numel(str2);
endInd = startInd+20;
RampupTime = str2num(header(startInd:endInd));

str = '<ParamLong."RampdownTime">';
ind = strfind(header, str);
len = numel(str);
str2 = '<Visible> "true"';
ind2 = strfind(header(ind+len:end), str2);
startInd = ind+len+ind2(1)+numel(str2);
endInd = startInd+20;
RampdownTime = str2num(header(startInd:endInd));

str = '<ParamLong."FlattopTime">';
ind = strfind(header, str);
len = numel(str);
str2 = '<Visible> "true"';
ind2 = strfind(header(ind+len:end), str2);
startInd = ind+len+ind2(1)+numel(str2);
endInd = startInd+20;
FlattopTime = str2num(header(startInd:endInd));

str = '<ParamLong."DelaySamplesTime">';
ind = strfind(header, str);
len = numel(str);
str2 = '<Visible> "true"';
ind2 = strfind(header(ind+len:end), str2);
startInd = ind+len+ind2(1)+numel(str2);
endInd = startInd+20;
DelaySamplesTime = str2num(header(startInd:endInd));

str = '<ParamDouble."ADCDuration">';
ind = strfind(header, str);
len = numel(str);
str2 = '<Precision> 16';
ind2 = strfind(header(ind+len:end), str2);
startInd = ind+len+ind2(1)+numel(str2);
endInd = startInd+20;
ADCDuration = str2num(header(startInd:endInd));

str = '<ParamLong."DestSamples">';
ind = strfind(header, str);
len = numel(str);
str2 = '<Visible> "true"';
ind2 = strfind(header(ind+len:end), str2);
startInd = ind+len+ind2(1)+numel(str2);
endInd = startInd+20;
DestSamples = str2num(header(startInd:endInd));
