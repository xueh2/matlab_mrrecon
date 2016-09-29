

function header = read_vb13_mdh(fid)

b                    = 0;

% typedef struct
% 
%  PACKED_MEMBER( uint32_t,     ulFlagsAndDMALength           );    // bit  0..24: DMA length [bytes]
%                                                                   // bit     25: pack bit
%                                                                   // bit 26..31: pci_rx enable flags                   4 byte
header.ulFlagsAndDMALength = fread(fid, 1,      'uint32');  b = b + 4;

%  PACKED_MEMBER( int32_t,      lMeasUID                      );    // measurement user ID                               4
header.lMeasUID            = fread(fid, 1,      'int32');   b = b + 4;

%  PACKED_MEMBER( uint32_t,     ulScanCounter                 );    // scan counter [1...]                               4
header.ulScanCounter       = fread(fid, 1,      'uint32');  b = b + 4;

%  PACKED_MEMBER( uint32_t,     ulTimeStamp                   );    // time stamp [2.5 ms ticks since 00:00]             4
header.ulTimeStamp         = fread(fid, 1,      'uint32');  b = b + 4;

%  PACKED_MEMBER( uint32_t,     ulPMUTimeStamp                );    // PMU time stamp [2.5 ms ticks since last trigger]  4
header.ulPMUTimeStamp      = fread(fid, 1,      'uint32');  b = b + 4;

%  PACKED_MEMBER( uint32_t,     aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK] ); // evaluation info mask field           8
header.aulEvalInfoMask     = fread(fid, [1 2],  'uint32');  b = b + 8;

%  PACKED_MEMBER( uint16_t,     ushSamplesInScan              );    // # of samples acquired in scan                     2
header.ushSamplesInScan    = fread(fid, 1,      'uint16'); b = b + 2;

%  PACKED_MEMBER( uint16_t,     ushUsedChannels               );    // # of channels used in scan                        2   =32
header.ushUsedChannels     = fread(fid, 1,      'uint16'); b = b + 2;

%  PACKED_MEMBER( sLoopCounter, sLC                           );    // loop counters                                    28   =60
%typedef struct
%
%  PACKED_MEMBER( uint16_t,  ushLine         ); /* line index                   */
%  PACKED_MEMBER( uint16_t,  ushAcquisition  ); /* acquisition index            */
%  PACKED_MEMBER( uint16_t,  ushSlice        ); /* slice index                  */
%  PACKED_MEMBER( uint16_t,  ushPartition    ); /* partition index              */
%  PACKED_MEMBER( uint16_t,  ushEcho         ); /* echo index                   */
%  PACKED_MEMBER( uint16_t,  ushPhase        ); /* phase index                  */
%  PACKED_MEMBER( uint16_t,  ushRepetition   ); /* measurement repeat index     */
%  PACKED_MEMBER( uint16_t,  ushSet          ); /* set index                    */
%  PACKED_MEMBER( uint16_t,  ushSeg          ); /* segment index  (for TSE)     */
%  PACKED_MEMBER( uint16_t,  ushIda          ); /* IceDimension a index         */
%  PACKED_MEMBER( uint16_t,  ushIdb          ); /* IceDimension b index         */
%  PACKED_MEMBER( uint16_t,  ushIdc          ); /* IceDimension c index         */
%  PACKED_MEMBER( uint16_t,  ushIdd          ); /* IceDimension d index         */
%  PACKED_MEMBER( uint16_t,  ushIde          ); /* IceDimension e index         */
% sLoopCounter;                                /* sizeof : 28 byte             */
header.sLC                        = fread(fid, [1 14], 'uint16'); b = b + 28;

%  PACKED_MEMBER( sCutOffData,  sCutOff                       );    // cut-off values                                    4
%typedef struct
%
%  PACKED_MEMBER( uint16_t,  ushPre          );    /* write ushPre zeros at line start */
%  PACKED_MEMBER( uint16_t,  ushPost         );    /* write ushPost zeros at line end  */
% sCutOffData;
header.sCutOff                    = fread(fid, [1 2],  'uint16'); b = b + 4;

%  PACKED_MEMBER( uint16_t,     ushKSpaceCentreColumn         );    // centre of echo                                    2
header.ushKSpaceCentreColumn      = fread(fid, 1,      'uint16'); b = b + 2;

%  PACKED_MEMBER( uint16_t,     ushCoilSelect                 );    // Bit 0..3: CoilSelect                              2
header.ushCoilSelect              = fread(fid, 1,      'uint16'); b = b + 2;

%  PACKED_MEMBER( float,        fReadOutOffcentre             );    // ReadOut offcenter value                           4
header.fReadOutOffcentre          = fread(fid, 1,      'float32');  b = b + 4;

%  PACKED_MEMBER( uint32_t,     ulTimeSinceLastRF             );    // Sequence time stamp since last RF pulse           4
header.ulTimeSinceLastRF          = fread(fid, 1,      'uint32');  b = b + 4;

%  PACKED_MEMBER( uint16_t,     ushKSpaceCentreLineNo         );    // number of K-space centre line                     2
header.ushKSpaceCentreLineNo      = fread(fid, 1,      'uint16'); b = b + 2;

%  PACKED_MEMBER( uint16_t,     ushKSpaceCentrePartitionNo    );    // number of K-space centre partition                2
header.ushKSpaceCentrePartitionNo = fread(fid, 1,      'uint16'); b = b + 2;

%  PACKED_MEMBER( uint16_t,     aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA] ); // free parameter for IceProgram   8   =88
header.aushIceProgramPara         = fread(fid, [1 4],  'uint16'); b = b + 8;

%  PACKED_MEMBER( uint16_t,     aushFreePara[MDH_FREEHDRPARA] );    // free parameter                          4 * 2 =   8
header.aushFreePara               = fread(fid, [1 4],  'uint16'); b = b + 8;

%  PACKED_MEMBER( sSliceData,   sSD                           );    // Slice Data                                       28   =124
%typedef struct
%
%  PACKED_MEMBER( float,  flSag          );
%  PACKED_MEMBER( float,  flCor          );
%  PACKED_MEMBER( float,  flTra          );
% sVector;
%typedef struct
%
%  PACKED_MEMBER( sVector,         sSlicePosVec     ); /* slice position vector        */
%  PACKED_MEMBER( float,           aflQuaternion[4] ); /* rotation matrix as quaternion*/
% sSliceData;                                         /* sizeof : 28 byte    */

header.sSD                        = fread(fid, [1 7],  'float32');  b = b + 28;

%  PACKED_MEMBER( uint16_t,	   ushPTABPosNeg                 );    // negative, absolute PTAB position in [0.1 mm]      2

%  PACKED_MEMBER( uint16_t,	   ushChannelId                  );    // channel Id must be the last parameter             2
header.ushChannelId               = fread(fid, 1,      'uint16');  b = b + 2;

%  sMDH;                                        // total length: 32 * 32 Bit (128 Byte)            128
% 
% #endif   /* MDH_H */

dummy                         = fread(fid, [1 128-b], 'uchar');
