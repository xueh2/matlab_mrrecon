function [COL, LIN, CHA, ACQ, SLC, PAR, ECO, PHS, REP, SET, SEG] = findKSpaceDimension(kspace)
% [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]

s = size(kspace);

COL = s(1);
LIN = s(2);
CHA = s(3);

ACQ = 1;
SLC = 1;
PAR = 1;
ECO = 1;
PHS = 1;
REP = 1;
SET = 1;
SEG = 1;

maxDimKSpace = 3;

if ( numel(s) >= 4 ) 
    ACQ = s(4); 
end

if ( numel(s) >= 5 ) 
    SLC = s(5); 
end

if ( numel(s) >= 6 ) PAR = s(6); end
if ( numel(s) >= 7 ) ECO = s(7); end
if ( numel(s) >= 8 ) PHS = s(8); end
if ( numel(s) >= 9 ) REP = s(9); end
if ( numel(s) >= 10 ) SET = s(10); end
if ( numel(s) >= 11 ) SEG = s(11); end
