function [header] = read_n4_header_ver4(datfilename, readoutindicator, readoutnumber);
% function [header] = read_n4_header_ver4(rawfilename, readoutindicator, readoutnumber);
%
% function to read siemens raw data file header for specific readoutnumber
%
% input:
%     datfilename   name of rawfile (version VB15/17 format)
%     readoutindicator   readoutindicator
%     readoutnumber number of readout specified for mdh header to output
% output:
%     header        matlab structure variable containing MDH header fields.

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************


position=readoutindicator(readoutnumber,12)-128;

fid=fopen(datfilename,'r','n');
fseek(fid,position,0);

header.ulDMALength         = fread(fid,1,'uint32'); % bit  0..27: DMA length [bytes]
                                                    % bit 28..31: pci_rx enable flags
mask=base2dec('0000111111111111111111111111111',2);
header.ulDMALength=bitand(header.ulDMALength,mask);
header.lMeasUID            = fread(fid,1,'int32'); % measurement user ID
header.ulScanCounter       = fread(fid,1,'uint32');  % scan counter [1...]
header.ulTimeStamp         = fread(fid,1,'uint32');  % time stamp [2.5 ms ticks since 00:00]
header.ulPMUTimeStamp      = fread(fid,1,'uint32'); % PMU time stamp [2.5 ms ticks since last trigger]
header.ulEvalInfoMask      = fread(fid,2,'uint32');  % evaluation info mask
header.ushSamplesInScan    = fread(fid,1,'uint16');  % # of samples acquired in scan
header.ushUsedChannels     = fread(fid,1,'uint16'); % # of channels used in scan

% loop counters
header.sLC.ushLine         = 1+fread(fid,1,'uint16'); % line index                   */
header.sLC.ushAcquisition  = 1+fread(fid,1,'uint16'); % acquisition index            */
header.sLC.ushSlice        = 1+fread(fid,1,'uint16'); % slice index                  */
header.sLC.ushPartition    = 1+fread(fid,1,'uint16'); % partition index              */
header.sLC.ushEcho         = 1+fread(fid,1,'uint16'); % echo index                   */
header.sLC.ushPhase        = 1+fread(fid,1,'uint16'); % phase index                  */
header.sLC.ushRepetition   = 1+fread(fid,1,'uint16'); % measurement repeat index     */
header.sLC.ushSet          = 1+fread(fid,1,'uint16'); % set index                    */
header.sLC.ushSeg          = 1+fread(fid,1,'uint16'); % segment index  (for TSE)     */
header.sLC.ushIda =fread(fid,1,'ushort');  % IceDimension a index        */
header.sLC.ushIdb =fread(fid,1,'ushort');  % IceDimension b index        */
header.sLC.ushIdc =fread(fid,1,'ushort');  % IceDimension c index        */
header.sLC.ushIdd =fread(fid,1,'ushort');  % IceDimension d index        */
header.sLC.ushIde =fread(fid,1,'ushort');  % IceDimension e index        */
% cut-off values  
header.sCutOffData.ushPre  = fread(fid,1,'uint16'); % write ushPre zeros at line start */
header.sCutOffData.ushPost = fread(fid,1,'uint16'); % write ushPost zeros at line end  */
header.ushKSpaceCentreColumn   = fread(fid,1,'uint16');  % centre of echo
header.ushCoilSelect       = fread(fid,1,'uint16');  % Bit 0..3: CoilSelect
header.fReadOutOffcentre       = fread(fid,1,'single');  % ReadOut offcenter value
header.ulTimeSinceLastRF       = fread(fid,1,'int32');   % Sequence time stamp since last RF pulse           
header.ushKSpaceCentreLineNo   = fread(fid,1,'uint16');  % number of K-space centre line
header.ushKSpaceCentrePartitionNo = fread(fid,1,'uint16');  % number of K-space centre partition
header.aushIceProgramPara=fread(fid,4,'uint16');
header.aushFreePara=fread(fid,4,'uint16');
% slice data (28 bytes)
header.sSD.sVector.flSag =  fread(fid,1,'single');%  slice position vector        */
header.sSD.sVector.flCor =  fread(fid,1,'single');%  slice position vector        */
header.sSD.sVector.flTra =  fread(fid,1,'single');%  slice position vector        */
header.sSD.aflQuaternion.a =  fread(fid,1,'single'); % rotation matrix as quaternion*/
header.sSD.aflQuaternion.b =  fread(fid,1,'single'); % rotation matrix as quaternion*/
header.sSD.aflQuaternion.c =  fread(fid,1,'single'); % rotation matrix as quaternion*/
header.sSD.aflQuaternion.d =  fread(fid,1,'single'); % rotation matrix as quaternion*/
header.ulChannelId = fread(fid,1,'int32');  % channel Id must be the last parameter

fclose(fid);

return
