
function [kspace, kspaceASCIndex, Noise, ref, refASCIndex, phsCorr, phsCorrASCIndex, reflect, reflectRef, reflectPhsCorr, other, otherASCIndex, SamplesInScan, KSpaceCentreColumn, MaxKSpaceLineNo, KSpaceCentreLineNo, KSpaceCentrePartitionNo] = Parse_ICE_Rawdata_VD(Data, asc) 
% read the Ice meas data
% /*--------------------------------------------------------------------------*/
% /* Definition of loop counter structure                                     */
% /* Note: any changes of this structure affect the corresponding swapping    */
% /*       method of the measurement data header proxy class (MdhProxy)       */
% /*--------------------------------------------------------------------------*/
% //-----------------------------------------------------------------------------
% for the VD line
% /*---------------------------------------------------------------------------*/
% /*  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential        */
% /*---------------------------------------------------------------------------*/
% /*
%  * Project: NUMARIS/4
%  *    File: \n4_servers1\pkg\MrServers\MrMeasSrv\SeqIF\MDH\mdh.h
%  * Version:
%  *  Author: CC_MEAS SCHOSTZF
%  *    Date: n.a.
%  *
%  *    Lang: C
%  *
%  * Descrip: measurement data header
%  *
%  *---------------------------------------------------------------------------*/
% 
% /*--------------------------------------------------------------------------*/
% /* Include control                                                          */
% /*--------------------------------------------------------------------------*/
% #ifndef MDH_H
% #define MDH_H
% 
% /*--------------------------------------------------------------------------*/
% /* Include MR basic type definitions                                        */
% /*--------------------------------------------------------------------------*/
% #include "MrCommon/MrGlobalDefinitions/MrBasicTypes.h"
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of header parameters                                         */
% /*--------------------------------------------------------------------------*/
% #define MDH_NUMBEROFEVALINFOMASK   2
% #define MDH_NUMBEROFICEPROGRAMPARA 24
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of free header parameters (short)                            */
% /*--------------------------------------------------------------------------*/
% #define MDH_RESERVEDHDRPARA  (4)
% 
% /*--------------------------------------------------------------------------*/
% /* Definition of time stamp tick interval/frequency                         */
% /* (used for ulTimeStamp and ulPMUTimeStamp                                 */
% /*--------------------------------------------------------------------------*/
% #define RXU_TIMER_INTERVAL  (2500000)     /* data header timer interval [ns]*/
% #define RXU_TIMER_FREQUENCY (400)         /* data header timer frequency[Hz]*/
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of bit masks for ulFlagsAndDMALength field                   */
% /*--------------------------------------------------------------------------*/
% #define MDH_DMA_LENGTH_MASK   (0x01FFFFFFL)
% #define MDH_PACK_BIT_MASK     (0x02000000L)
% #define MDH_ENABLE_FLAGS_MASK (0xFC000000L)
% 
% /*--------------------------------------------------------------------------*/
% /* Definition of loop counter structure                                     */
% /* Note: any changes of this structure affect the corresponding swapping    */
% /*       method of the measurement data header proxy class (MdhProxy)       */
% /*--------------------------------------------------------------------------*/
% #include "MrServers/MrVista/include/pack.h"
% 
% /// \ingroup MDH
% /// \todo write documentation 
% /// \brief Definition of loop counter structure
% typedef struct
% {
%   PACKED_MEMBER( uint16_t,  ushLine         ); /**< line index                   */
%   PACKED_MEMBER( uint16_t,  ushAcquisition  ); /**< acquisition index            */
%   PACKED_MEMBER( uint16_t,  ushSlice        ); /**< slice index                  */
%   PACKED_MEMBER( uint16_t,  ushPartition    ); /**< partition index              */
%   PACKED_MEMBER( uint16_t,  ushEcho         ); /**< echo index                   */
%   PACKED_MEMBER( uint16_t,  ushPhase        ); /**< phase index                  */
%   PACKED_MEMBER( uint16_t,  ushRepetition   ); /**< measurement repeat index     */
%   PACKED_MEMBER( uint16_t,  ushSet          ); /**< set index                    */
%   PACKED_MEMBER( uint16_t,  ushSeg          ); /**< segment index  (for TSE)     */
%   PACKED_MEMBER( uint16_t,  ushIda          ); /**< IceDimension a index         */
%   PACKED_MEMBER( uint16_t,  ushIdb          ); /**< IceDimension b index         */
%   PACKED_MEMBER( uint16_t,  ushIdc          ); /**< IceDimension c index         */
%   PACKED_MEMBER( uint16_t,  ushIdd          ); /**< IceDimension d index         */
%   PACKED_MEMBER( uint16_t,  ushIde          ); /**< IceDimension e index         */
% } sLoopCounter;                                /* sizeof : 28 byte             */
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of slice vectors                                             */
% /*--------------------------------------------------------------------------*/
% 
% /// \ingroup MDH
% /// \todo write documentation 
% /// \brief Definition of slice vectors 
% typedef struct
% {
%   PACKED_MEMBER( float,  flSag          );
%   PACKED_MEMBER( float,  flCor          );
%   PACKED_MEMBER( float,  flTra          );
% } sVector; /* 12 bytes */
% 
% /// \ingroup MDH
% /// \todo write documentation 
% /// \brief Definition of slice data structure
% typedef struct
% {
%   PACKED_STRUCT( sVector,         sSlicePosVec     ); /**< slice position vector        */
%   PACKED_MEMBER( float,           aflQuaternion[4] ); /**< rotation matrix as quaternion*/
% } sSliceData;                                         /* sizeof : 28 byte             */
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of cut-off data                                              */
% /*--------------------------------------------------------------------------*/
% /// \ingroup MDH
% /// \todo write documentation 
% /// \brief Definition of cut-off data
% typedef struct
% {
%   PACKED_MEMBER( uint16_t,  ushPre          );    /**< write ushPre zeros at line start */
%   PACKED_MEMBER( uint16_t,  ushPost         );    /**< write ushPost zeros at line end  */
% } sCutOffData; /* 4 bytes */
% 
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of measurement data header                                   */
% /*--------------------------------------------------------------------------*/
% /// \ingroup MDH
% /// \todo write documentation 
% /// \brief Definition of the scan header structure
% typedef struct sScanHeader
% {
%   PACKED_MEMBER( uint32_t,     ulFlagsAndDMALength           );                 ///<  0: ( 4) bit  0..24: DMA length [bytes]
%                                                                                 ///<          bit     25: pack bit
%                                                                                 ///<          bit 26..31: pci_rx enable flags
%   PACKED_MEMBER( int32_t,      lMeasUID                      );                 ///<  4: ( 4) measurement user ID
%   PACKED_MEMBER( uint32_t,     ulScanCounter                 );                 ///<  8: ( 4) scan counter [1...]
%   PACKED_MEMBER( uint32_t,     ulTimeStamp                   );                 ///< 12: ( 4) time stamp [2.5 ms ticks since 00:00]
%   PACKED_MEMBER( uint32_t,     ulPMUTimeStamp                );                 ///< 16: ( 4) PMU time stamp [2.5 ms ticks since last trigger]
%   PACKED_MEMBER( uint16_t,     ushSystemType                 );                 ///< 20: ( 2) System type (todo: values?? ####)
%   PACKED_MEMBER( uint16_t,     ulPTABPosDelay                );                 ///< 22: ( 2) PTAb delay ??? TODO: How do we handle this ####
%   PACKED_MEMBER( int32_t,	     lPTABPosX                     );                 ///< 24: ( 4) absolute PTAB position in [µm]
%   PACKED_MEMBER( int32_t,	     lPTABPosY                     );                 ///< 28: ( 4) absolute PTAB position in [µm]
%   PACKED_MEMBER( int32_t,	     lPTABPosZ                     );                 ///< 32: ( 4) absolute PTAB position in [µm]
%   PACKED_MEMBER( uint32_t,	   ulReserved1                   );                 ///< 36: ( 4) reserved for future hardware signals
%   PACKED_MEMBER( uint32_t,     aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK]);      ///< 40: ( 8) evaluation info mask field
%   PACKED_MEMBER( uint16_t,     ushSamplesInScan              );                 ///< 48: ( 2) # of samples acquired in scan
%   PACKED_MEMBER( uint16_t,     ushUsedChannels               );                 ///< 50: ( 2) # of channels used in scan
%   PACKED_STRUCT( sLoopCounter, sLC                           );                 ///< 52: (28) loop counters
%   PACKED_STRUCT( sCutOffData,  sCutOff                       );                 ///< 80: ( 4) cut-off values
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentreColumn         );                 ///< 84: ( 2) centre of echo
%   PACKED_MEMBER( uint16_t,     ushCoilSelect                 );                 ///< 86: ( 2) Bit 0..3: CoilSelect
%   PACKED_MEMBER( float,        fReadOutOffcentre             );                 ///< 88: ( 4) ReadOut offcenter value
%   PACKED_MEMBER( uint32_t,     ulTimeSinceLastRF             );                 ///< 92: ( 4) Sequence time stamp since last RF pulse
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentreLineNo         );                 ///< 96: ( 2) number of K-space centre line
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentrePartitionNo    );                 ///< 98: ( 2) number of K-space centre partition
%   PACKED_STRUCT( sSliceData,   sSD                           );                 ///< 100:(28) Slice Data
%   PACKED_MEMBER( uint16_t,     aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA] );///< 128:(48) free parameter for IceProgram
%   PACKED_MEMBER( uint16_t,     aushReservedPara[MDH_RESERVEDHDRPARA] );         ///< 176:( 8) unused parameter (padding to next 192byte alignment )
%                                                                                 ///<          NOTE: These parameters MUST NOT be used by any application (for future use)
%   PACKED_MEMBER( uint16_t,     ushApplicationCounter         );                 ///< 184 ( 2)
%   PACKED_MEMBER( uint16_t,     ushApplicationMask            );                 ///< 186 ( 2)
%   PACKED_MEMBER( uint32_t,     ulCRC                         );                 ///< 188:( 4) CRC 32 checksum
% } sScanHeader;                                                                  // total length: 6 x 32 Byte (192 Byte)
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of channel data header                                   */
% /*--------------------------------------------------------------------------*/
% /// \ingroup MDH
% /// \todo write documentation 
% /// \brief Definition of the scan header structure
% typedef struct sChannelHeader
% {
%   PACKED_MEMBER( uint32_t,     ulTypeAndChannelLength        );    ///< 0: (4) bit  0.. 7: type (0x02 => ChannelHeader)
%                                                                    ///<        bit  8..31: channel length (header+data) in byte
%                                                                    ///<        type   := ulTypeAndChannelLength & 0x000000FF
%                                                                    ///<        length := ulTypeAndChannelLength >> 8
%   PACKED_MEMBER( int32_t,      lMeasUID                      );    ///< 4: (4) measurement user ID
%   PACKED_MEMBER( uint32_t,     ulScanCounter                 );    ///< 8: (4) scan counter [1...]
%   PACKED_MEMBER( uint32_t,     ulReserved1                   );    ///< 12:(4) reserved
%   PACKED_MEMBER( uint32_t,     ulSequenceTime                );    ///< 16:(4) Sequence readout starting time bit 31..9 time in [10us]
%                                                                    ///<                                       bit  8..0 time in [25ns]
%   PACKED_MEMBER( uint32_t,     ulUnused2                     );    ///< 20:(4) unused
%   PACKED_MEMBER( uint16_t,     ulChannelId                   );    ///< 24:(4) unused
%   PACKED_MEMBER( uint16_t,     ulUnused3                     );    ///< 26:(2) unused
%   PACKED_MEMBER( uint32_t,     ulCRC                         );    ///< 28:(4) CRC32 checksum of channel header
% } sChannelHeader;                                                  // total length:  32 byte
% 
% #include "MrServers/MrVista/include/unpack.h"
% 
% #endif   /* MDH_H */
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of EvalInfoMask:                                             */
% /*--------------------------------------------------------------------------*/
% /// \ingroup MDH
% /// \todo write documentation 
% /// \brief Bits in the EvalInfoMask
% /// \see sMDH
% /// @{
% const MdhBitField MDH_ACQEND            (0UL);       ///< last scan 
% const MdhBitField MDH_RTFEEDBACK        (1UL);       ///< Realtime feedback scan
% const MdhBitField MDH_HPFEEDBACK        (2UL);       ///< High perfomance feedback scan
% const MdhBitField MDH_ONLINE            (3UL);       ///< processing should be done online
% const MdhBitField MDH_OFFLINE           (4UL);       ///< processing should be done offline
% const MdhBitField MDH_SYNCDATA          (5UL);       ///< readout contains synchroneous data
% const MdhBitField MDH_LASTSCANINCONCAT  (8UL);       ///< Flag for last scan in concatination
% 
% const MdhBitField MDH_RAWDATACORRECTION (10UL);      ///< Correct the rawadata with the rawdata correction factor
% const MdhBitField MDH_LASTSCANINMEAS    (11UL);      ///< Flag for last scan in measurement
% const MdhBitField MDH_SCANSCALEFACTOR   (12UL);      ///< Flag for scan specific additional scale factor
% const MdhBitField MDH_2NDHADAMARPULSE   (13UL);      ///< 2nd RF exitation of HADAMAR
% const MdhBitField MDH_REFPHASESTABSCAN  (14UL);      ///< reference phase stabilization scan
% const MdhBitField MDH_PHASESTABSCAN     (15UL);      ///< phase stabilization scan
% const MdhBitField MDH_D3FFT             (16UL);      ///< execute 3D FFT
% const MdhBitField MDH_SIGNREV           (17UL);      ///< sign reversal
% const MdhBitField MDH_PHASEFFT          (18UL);      ///< execute phase fft
% const MdhBitField MDH_SWAPPED           (19UL);      ///< swapped phase/readout direction
% const MdhBitField MDH_POSTSHAREDLINE    (20UL);      ///< shared line
% const MdhBitField MDH_PHASCOR           (21UL);      ///< phase correction data
% const MdhBitField MDH_PATREFSCAN        (22UL);      ///< additonal scan for PAT reference line/partition
% const MdhBitField MDH_PATREFANDIMASCAN  (23UL);      ///< additonal scan for PAT reference line/partition that is also used as image scan
% const MdhBitField MDH_REFLECT           (24UL);      ///< reflect line
% const MdhBitField MDH_NOISEADJSCAN      (25UL);      ///< noise adjust scan 
% const MdhBitField MDH_SHARENOW          (26UL);      ///< all lines are acquired from the actual and previous e.g. phases
% const MdhBitField MDH_LASTMEASUREDLINE  (27UL);      ///< indicates that the current line is the last measured line of all succeeding e.g. phases
% const MdhBitField MDH_FIRSTSCANINSLICE  (28UL);      ///< indicates first scan in slice (needed for time stamps)
% const MdhBitField MDH_LASTSCANINSLICE   (29UL);      ///< indicates  last scan in slice (needed for time stamps)
% const MdhBitField MDH_TREFFECTIVEBEGIN  (30UL);      ///< indicates the begin time stamp for TReff (triggered measurement)
% const MdhBitField MDH_TREFFECTIVEEND    (31UL);      ///< indicates the   end time stamp for TReff (triggered measurement)
% const MdhBitField MDH_MDS_REF_POSITION  (32UL);      ///< indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
% const MdhBitField MDH_SLC_AVERAGED      (33UL);      ///< indicates avveraged slice for slice partial averaging scheme
% const MdhBitField MDH_TAGFLAG1          (34UL);      ///< adjust scan 
% 
% const MdhBitField MDH_CT_NORMALIZE              (35UL);  ///< Marks scans used to calculate correction maps for TimCT-Prescan normalize
% const MdhBitField MDH_SCAN_FIRST                (36UL);  ///< Marks the first scan of a particular map
% const MdhBitField MDH_SCAN_LAST                 (37UL);  ///< Marks the last scan of a particular map
% 
% const MdhBitField MDH_FIRST_SCAN_IN_BLADE       (40UL);  ///< Marks the first line of a blade
% const MdhBitField MDH_LAST_SCAN_IN_BLADE        (41UL);  ///< Marks the last line of a blade
% const MdhBitField MDH_LAST_BLADE_IN_TR          (42UL);  ///< Set for all lines of the last BLADE in each TR interval
%                                                           
% const MdhBitField MDH_PACE                      (44UL);  ///< Distinguishes PACE scans from non PACE scans.
%                                                           
% const MdhBitField MDH_RETRO_LASTPHASE           (45UL);  ///< Marks the last phase in a heartbeat
% const MdhBitField MDH_RETRO_ENDOFMEAS           (46UL);  ///< Marks an ADC at the end of the measurement
% const MdhBitField MDH_RETRO_REPEATTHISHEARTBEAT (47UL);  ///< Repeat the current heartbeat when this bit is found
% const MdhBitField MDH_RETRO_REPEATPREVHEARTBEAT (48UL);  ///< Repeat the previous heartbeat when this bit is found
% const MdhBitField MDH_RETRO_ABORTSCANNOW        (49UL);  ///< Just abort everything
% const MdhBitField MDH_RETRO_LASTHEARTBEAT       (50UL);  ///< This adc is from the last heartbeat (a dummy)
% const MdhBitField MDH_RETRO_DUMMYSCAN           (51UL);  ///< This adc is just a dummy scan, throw it away
% const MdhBitField MDH_RETRO_ARRDETDISABLED      (52UL);  ///< Disable all arrhythmia detection when this bit is found
% /// @}
% 
% //-----------------------------------------------------------------------------
% // Definition of EvalInfoMask for COP:
% //-----------------------------------------------------------------------------
% #define	MDH_COP_ACQEND            (0x00000001)
% #define MDH_COP_RTFEEDBACK        (0x00000002)
% #define MDH_COP_HPFEEDBACK        (0x00000004)
% #define	MDH_COP_ONLINE            (0x00000008)
% #define	MDH_COP_OFFLINE           (0x00000010)
% #define	MDH_COP_SYNCDATA          (0x00000020)
% /*---------------------------------------------------------------------------*/
% /*  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential        */
% /*---------------------------------------------------------------------------*/
% //-----------------------------------------------------------------------------
% [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
% 11 dimension array
%
% 
%     enum ReorderingScheme
%     {
%         RS_LIN_IN_PAR          = 1,   // lines in partition
%         RS_PAR_IN_LIN          = 2,   // partitions in lines
%         RS_ARBITRARY           = 3    // arbitrary reordering schemes e.g. centric reordering
%     };
%     
%     enum PartialFourierFactor
%     {
%         PF_HALF = 0x01,
%         PF_5_8  = 0x02,
%         PF_6_8  = 0x04,
%         PF_7_8  = 0x08,
%         PF_OFF  = 0x10,
%         PF_AUTO = 0x20
%     };
%     
%     enum AveragingMode
%     {
%         INNER_LOOP = 0x01,
%         OUTER_LOOP = 0x02
%     };
%     
%     enum MultiSliceMode
%     {
%         MSM_SEQUENTIAL  = 0x01,
%         MSM_INTERLEAVED = 0x02,
%         MSM_SINGLESHOT  = 0x04
%     };
%     
%     enum Dimension
%     {
%         DIM_1 = 0x01,
%         DIM_2 = 0x02,
%         DIM_3 = 0x04
%     };
%     
%     enum Trajectory
%     {
%         TRAJECTORY_CARTESIAN  = 0x01,
%         TRAJECTORY_RADIAL     = 0x02,
%         TRAJECTORY_SPIRAL     = 0x04,
%         TRAJECTORY_BLADE      = 0x08,
%     };
%     
%     enum ViewSharing
%     {
%         VIEW_SHARING_OFF       = 0x01,
%         VIEW_SHARING_SHPHS     = 0x02,
%         VIEW_SHARING_SHVENC    = 0x04,
%         VIEW_SHARING_TWIST     = 0x08
%     };
%     
%     enum AsymmetricEchoMode
%     {
%         ASYMM_ECHO_WEAK       = 0x1,
%         ASYMM_ECHO_STRONG     = 0x2,
%         ASYMM_ECHO_HALF       = 0x4
%     };
% 
%     enum Reordering
%     {
%         REORDERING_LINEAR         = 0x01,
%         REORDERING_CENTRIC        = 0x02,
%         REORDERING_LINE_SEGM      = 0x04,
%         REORDERING_PART_SEGM      = 0x08,
%         REORDERING_FREE_0         = 0x10,
%         REORDERING_FREE_1         = 0x20,
%         REORDERING_FREE_2         = 0x40,
%         REORDERING_FREE_3         = 0x80,
%         REORDERING_LINEAR_ROTATED = 0x0100,
%         REORDERING_RADIAL         = 0x0200,
%         REORDERING_SPIRAL         = 0x0400       
%     };
% 
%     enum POCS //(projection onto convex sets recon mode)
%     {
%         POCS_OFF        = 0x01,
%         POCS_READ_SLICE = 0x02,
%         POCS_READ_PHASE = 0x04
%     };
%     
%     enum Increment
%     {
%         INC_NORMAL           = 1,
%         INC_BASE2            = 2,
%         INC_FIX              = 3,
%         INC_TURBO_FACTOR     = 4,
%         INC_EPI_FACTOR       = 5,
%         INC_SEGMENTED        = 6,
%         INC_TGSE_FACTOR      = 7,
%         INC_GRE_SEGMENTS     = 8,
%         INC_SINGLESHOT       = 9,
%         INC_TWICE_EPI_FACTOR = 10,
%         INC_64               = 11,
%         INC_32               = 12,
%         INC_16               = 13,
%         /// UILink / KSpace taking care of Reflines and IPAT
%         INC_NORMAL_IPAT           = 14,
%         INC_BASE2_IPAT            = 15,
%         INC_FIX_IPAT              = 16,
%         INC_TURBO_FACTOR_IPAT     = 17,
%         INC_EPI_FACTOR_IPAT       = 18,
%         INC_SEGMENTED_IPAT        = 19,
%         INC_TGSE_FACTOR_IPAT      = 20,
%         INC_GRE_SEGMENTS_IPAT     = 21,
%         INC_SINGLESHOT_IPAT       = 22,
%         INC_TWICE_EPI_FACTOR_IPAT = 23,
%         INC_64_IPAT               = 24,
%         INC_32_IPAT               = 25,
%         INC_16_IPAT               = 26,
%     };
% };
% 
% #endif


tic
S0 = size(Data);
if ( numel(S0) > 2 )
    Data = reshape(Data, [S0(1) S0(2)*S0(3)]);
end

S0 = size(Data);
maxCOL = S0(1);
numOfLines = S0(2);
numOfCHA = double(asc(1).ushUsedChannels);
numOfScans = numOfLines/numOfCHA;

maxCOL = 0;
startInd = double(asc(1).ushUsedChannels+1);
for ii=startInd:numOfLines    
    [bNoiseLine, bSeperateRef, bRefandImg, bPhaseCorr, bReflect, bLastScanInSlice, bRTFeedBack] = parseEvalInfoMask(asc(ii).aulEvalInfoMask(1));
    if ( ~bNoiseLine && ~bRTFeedBack && asc(ii).ushSamplesInScan > maxCOL )
        maxCOL = asc(ii).ushSamplesInScan;
    end
end

% find the noise lines and seperate ref lines if any
numOfNoiseLines = 0; 
NoiseLines = zeros(numOfScans,1, 'int32');
sLCNoise = zeros(numOfScans, 11, 'int32');

numOfSeperateRefLines = 0; 
RefLines = zeros(numOfScans,1, 'int32');
numOfRefAndImgLines = 0; 
RefAndImgLines = zeros(numOfScans,1, 'int32');
sLCRef = zeros(numOfScans, 11, 'int32');

numOfReflectRefLines = 0; 
ReflectRefLines = zeros(numOfScans,1, 'int32');
sLCReflectRef = zeros(numOfScans, 11, 'int32');

numOfPhaseCorrLines = 0; 
PhaseCorrLines = zeros(numOfScans,1, 'int32');
sLCPhaseCorr = zeros(numOfScans, 11, 'int32');

numOfReflectLines = 0; 
ReflectLines = zeros(numOfScans,1, 'int32');
sLCReflect = zeros(numOfScans, 11, 'int32');

numOfDataLines = 0; 
DataLines = zeros(numOfScans,1, 'int32');
sLC = zeros(numOfScans, 11, 'int32');

numOfOtherLines = 0; 
OtherLines = zeros(numOfScans,1, 'int32');
sLCOther = zeros(numOfScans, 11, 'int32');

kspaceCentreLineNo = zeros(numOfScans,1, 'int32');

currREPForOther = 1;
for s=1:numOfScans
    
    % starting line for this scan
    i = (s-1)*numOfCHA + 1;
    
    if ( IsAcqEnd(asc(i).aulEvalInfoMask(1)) )
        continue;
    end
    
    kspaceCentreLineNo(s) = asc(i).ushKSpaceCentreLineNo;
    
    [bNoiseLine, bSeperateRef, bRefandImg, bPhaseCorr, bReflect, bLastScanInSlice, bRTFeedBack] = parseEvalInfoMask(asc(i).aulEvalInfoMask(1));
    
    if ( bRTFeedBack )
        continue;
    end
    
    if ( bNoiseLine  )
        numOfNoiseLines = numOfNoiseLines + 1;
        NoiseLines(numOfNoiseLines) = s;
        sLCNoise(s,:) = fillSLC(asc(i));
        continue;
    end
    
    if ( bPhaseCorr )
        numOfPhaseCorrLines = numOfPhaseCorrLines + 1;
        PhaseCorrLines(numOfPhaseCorrLines) = s;
        sLCPhaseCorr(s,:) = fillSLC(asc(i));
        if ( bReflect ) % if it is a reflected phase correction lines
            numOfReflectLines = numOfReflectLines + 1;
            ReflectLines(numOfReflectLines) = s;
            sLCReflect(s,:) = fillSLC(asc(i));
        end
        continue;
    end 
    
    if ( bSeperateRef || bRefandImg )
        numOfSeperateRefLines = numOfSeperateRefLines + 1;
        RefLines(numOfSeperateRefLines) = s;
        sLCRef(s,:) = fillSLC(asc(i));
        
        if ( bReflect )
            numOfReflectRefLines = numOfReflectRefLines + 1;
            ReflectRefLines(numOfReflectRefLines) = s;
            sLCReflectRef(s,:) = fillSLC(asc(i));
        end
        
        if ( bSeperateRef & ~bRefandImg )
            continue;
        end
        
        if ( bRefandImg )
            numOfRefAndImgLines = numOfRefAndImgLines + 1;
            RefAndImgLines(numOfRefAndImgLines) = s;
        end
    end
        
    if ( asc(i).ushSamplesInScan < maxCOL & ~bRefandImg  )
        numOfOtherLines = numOfOtherLines + 1;
        OtherLines(numOfOtherLines) = s;
        sLCOther(s,:) = fillSLC(asc(i));
        sLCOther(s,9) = currREPForOther;
        continue;
    end
        
    if ( bReflect ) % if it is a reflected data lines
        numOfReflectLines = numOfReflectLines + 1;
        ReflectLines(numOfReflectLines) = s;
        sLCReflect(s,:) = fillSLC(asc(i));
    end 

    numOfDataLines = numOfDataLines + 1; % include normal data line and not refAndImg lines
    DataLines(numOfDataLines) = s;
    sLC(s,:) = fillSLC(asc(i));

    if ( bLastScanInSlice & sLC(numOfDataLines,5)==1 )
        currREPForOther = sLC(numOfDataLines,9) + 1;
        % disp(['current REP : ' num2str(currREPForOther)]);
    end
end

NoiseLines = NoiseLines(1:numOfNoiseLines);
RefLines = RefLines(1:numOfSeperateRefLines); 
RefAndImgLines = RefAndImgLines(1:numOfRefAndImgLines);
ReflectRefLines = ReflectRefLines(1:numOfReflectRefLines);
PhaseCorrLines = PhaseCorrLines(1:numOfPhaseCorrLines);
ReflectLines = ReflectLines(1:numOfReflectLines);
DataLines = DataLines(1:numOfDataLines);
OtherLines = OtherLines(1:numOfOtherLines);

% % combination along SEG, for the sake of EPI recon, we shall not combine the SEG
% sLCRef(:,11) = 1;
% sLCPhaseCorr(:,11) = 1;
% sLCReflect(:,11) = 1;
% sLC(:,11) = 1;
% sLCOther(:,11) = 1;

maxKSpaceLineNo = 2*max(kspaceCentreLineNo);
[kspace, kspaceASCIndex] = allocateKSpace(sLC, maxKSpaceLineNo);
kspace(:) = 0;

MaxKSpaceLineNo = max(sLC(:,2));

% SEG = size(kspace, 11);

% noise data if any
if ( numOfNoiseLines > 0 )
    maxColNoise = max(sLCNoise(:, 1));
    maxLinNoise = max(sLCNoise(:, 2));
    maxChaNoise = max(sLCNoise(:, 3));
    Noise = zeros(maxColNoise, maxLinNoise, maxChaNoise);
else
    Noise = [];
end

% seperate reference lines
if ( numOfSeperateRefLines > 0 )
    % [ref, refASCIndex] = allocateKSpace(sLCRef, maxKSpaceLineNo);
    [ref, refASCIndex] = allocateKSpace(sLC, maxKSpaceLineNo);
    ref(:) = 0;
%     sRef = size(ref);
%     NRef = prod(sRef(2:11));
%     ref = reshape(ref, [sRef(1) NRef]);    
else
    ref = [];
    refASCIndex = [];
end

% reflect ref lines
if ( numOfReflectRefLines > 0 )
    ss = size(ref);    
    reflectRef = zeros([1 ss(2:end)], 'int8');
else
    reflectRef = [];
end

% phase correction lines
if ( numOfPhaseCorrLines > 0 )
    [phsCorr, phsCorrASCIndex] = allocateKSpace(sLCPhaseCorr, maxKSpaceLineNo);
    ss = size(phsCorr);
    reflectPhsCorr = zeros([1 ss(2:end)], 'int8');
    hitCountPhsCorr = zeros([1 ss(2:end)], 'int8');
else
    phsCorr = [];
    phsCorrASCIndex = [];
    reflectPhsCorr = [];
    hitCountPhsCorr = [];
end

% phase correction lines
if ( numOfReflectLines > 0 )
    ss = size(kspace);
    reflect = zeros([1 ss(2:end)], 'int8');
else
    reflect = [];
end

% other lines
if ( numOfOtherLines > 0 )
    [other, otherASCIndex] = allocateKSpace(sLCOther, maxKSpaceLineNo);
else
    other = [];
    otherASCIndex = [];
end

% put the data into the array

SamplesInScan = 0;
KSpaceCentreColumn = 0;
KSpaceCentreLineNo = 0;
KSpaceCentrePartitionNo = 0;

scanN = numel(NoiseLines);
if ( scanN > 0 )    
    for s=1:scanN
        i = double((NoiseLines(s)-1)*numOfCHA+1);
        N = asc(i).ushSamplesInScan;
        Noise(1:N, asc(i).sLC(1)+1, 1:numOfCHA) = Data(1:N, i:i+numOfCHA-1);
    end    
end

scanN = numel(PhaseCorrLines);
if ( scanN > 0 )    
    for s=1:scanN
        r = PhaseCorrLines(s);
        i = double((r-1)*numOfCHA+1);
        N = asc(i).ushSamplesInScan;
        
        if ( s==1 )
            aZ = addPrePostZeros(asc(i));
        end
        
        if ( sLCPhaseCorr(r,2) <= maxKSpaceLineNo )
            if ( aZ == 0 )
                phsCorr(1:N, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i:i+numOfCHA-1);
            elseif ( aZ == 1 ) % pre zeros
                phsCorr(end-N+1:end, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i:i+numOfCHA-1);
            elseif ( aZ == 2 ) % post zeros
                phsCorr(1:N, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i:i+numOfCHA-1);
            end
            phsCorrASCIndex(sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = i:i+numOfCHA-1;
                
            hitCountPhsCorr(1, sLCPhaseCorr(r,2), 1, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = hitCountPhsCorr(1, sLCPhaseCorr(r,2), 1, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) + 1;
        end

        % check whether it is a reflected line
        r2 = find(r==ReflectLines);
        if ( ~isempty(r2)  )
             reflectPhsCorr(1, sLCPhaseCorr(r,2), 1:numOfCHA, sLCReflect(r,4), sLCReflect(r,5), ...
                    sLCReflect(r,6), sLCReflect(r,7), sLCReflect(r,8), sLCReflect(r,9), sLCReflect(r,10), sLCReflect(r,11)) = 1;
        end
    end
    
    for s=1:scanN
        r = PhaseCorrLines(s);
        hitCount = hitCountPhsCorr(1, sLCPhaseCorr(r,2), 1, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11));
        
        if ( hitCount > 0 )
            phsCorr(:, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                        sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = phsCorr(:, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                        sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) ./ double(hitCount);
        end
    end
end
    
scanN = numel(RefLines);
if ( scanN > 0 )    
    for s=1:scanN
        
        if ( mod(s, 100) == 0 )
            disp(['ref scan ' num2str(s)]);
        end
        
        r = RefLines(s);
        i = double((r-1)*numOfCHA+1);
        
        if ( s == 1 )
            N = asc(i).ushSamplesInScan;
            aZ = addPrePostZeros(asc(i));
            
            if ( aZ == 0 )
                sFe = 1;
                eFe = N;
            elseif ( aZ == 1 ) % pre zeros
                sFe = size(ref,1)-N+1;
                eFe = size(ref,1);
            elseif ( aZ == 2 ) % post zeros
                sFe = 1;
                eFe = N;
            end
        end
        
        if ( sLCRef(r,2) <= maxKSpaceLineNo )
                        
%             if ( aZ == 0 )
%                 ref(1:N, sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 1 ) % pre zeros
%                 ref(end-N+1:end, sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 2 ) % post zeros
%                 ref(1:N, sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             end
            
            ref(sFe:eFe, sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i:i+numOfCHA-1);
                
            refASCIndex(sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
                    sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = i:i+numOfCHA-1;
        end      
    end    
end

scanN = numel(ReflectRefLines);
if ( scanN > 0 )    
    for s=1:scanN
        % check whether it is a reflected line
        r2 = ReflectRefLines(s);
        if ( ~isempty(r2)  )
             reflectRef(1, sLCReflectRef(r2,2), 1:numOfCHA, sLCReflectRef(r2,4), sLCReflectRef(r2,5), ...
                    sLCReflectRef(r2,6), sLCReflectRef(r2,7), sLCReflectRef(r2,8), sLCReflectRef(r2,9), sLCReflectRef(r2,10), sLCReflectRef(r2,11)) = 1;
        end
    end
end

scanN = numel(DataLines);
if ( scanN > 0 )    
    for s=1:scanN
        
        if ( mod(s, 100) == 0 )
            disp(['data scan ' num2str(s)]);
        end
        
        r = DataLines(s);
        l = r;
        
        i = double((r-1)*numOfCHA+1);
        
        N = asc(i).ushSamplesInScan;
        
        if ( SamplesInScan < N ) SamplesInScan = N; end
        if ( KSpaceCentreColumn < asc(i).ushKSpaceCentreColumn ) KSpaceCentreColumn = asc(i).ushKSpaceCentreColumn; end
        if ( KSpaceCentreLineNo < asc(i).ushKSpaceCentreLineNo ) KSpaceCentreLineNo = asc(i).ushKSpaceCentreLineNo; end
        if ( KSpaceCentrePartitionNo < asc(i).ushKSpaceCentrePartitionNo ) KSpaceCentrePartitionNo = asc(i).ushKSpaceCentrePartitionNo; end
        
        if ( s == 1 )
            aZ = addPrePostZeros(asc(i));
        end
        
        if ( sLC(l,2) <= maxKSpaceLineNo )
            if ( aZ == 0 )
                kspace(1:N, sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i:i+numOfCHA-1);
            elseif ( aZ == 1 ) % pre zeros
                kspace(end-N+1:end, sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i:i+numOfCHA-1);
            elseif ( aZ == 2 ) % post zeros
                kspace(1:N, sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i:i+numOfCHA-1);
            end
            
            kspaceASCIndex(sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = i:i+numOfCHA-1;
        end       
    end
end
    
scanN = numel(ReflectLines);
if ( scanN > 0 )    
    for s=1:scanN
        % check whether it is a reflected line
        r2 = ReflectLines(s);
        if ( ~isempty(r2)  )
             reflect(1, sLCReflect(r2,2), 1:numOfCHA, sLCReflect(r2,4), sLCReflect(r2,5), ...
                    sLCReflect(r2,6), sLCReflect(r2,7), sLCReflect(r2,8), sLCReflect(r2,9), sLCReflect(r2,10), sLCReflect(r2,11)) = 1;
        end
    end
end
    
scanN = numel(OtherLines);
if ( scanN > 0 )    
    for s=1:scanN
        
        r = OtherLines(s);
        i = double((r-1)*numOfCHA+1);
        N = asc(i).ushSamplesInScan;
        
        if ( s==1 )
            aZ = addPrePostZeros(asc(i));
        end
        
        if ( sLCOther(r,2) <= maxKSpaceLineNo )
                       
            if ( aZ == 0 )
                other(1:N, sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i:i+numOfCHA-1);
            elseif ( aZ == 1 ) % pre zeros
                other(end-N+1:end, sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i:i+numOfCHA-1);
            elseif ( aZ == 2 ) % post zeros
                other(1:N, sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i:i+numOfCHA-1);
            end
            
            otherASCIndex(sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = i:i+numOfCHA-1;
        end
    end
end

% l = 0;
% r = 0;
% n = 0;
% for s=1:numOfScans
%     s;
%     i = double((s-1)*numOfCHA+1);
%     
%     if ( mod(s, 1e3) == 0 )
%         disp([num2str(s) ' scans parsed ... ']);
%     end
%     
%     if ( IsAcqEnd(asc(i).aulEvalInfoMask(1)) )
%         continue;
%     end
% 
%     % fill noise lines
% %     r = find(s==NoiseLines);
% %     if ( ~isempty(r) )
% %         n = n + 1;
% %         N = asc(i).ushSamplesInScan;
% %         aZ = addPrePostZeros(asc(i));
% %         aZ = 0;
% %         if ( aZ == 0 )
% %             Noise(1:N, sLCNoise(r,2), 1:numOfCHA) = Data(1:N, i:i+numOfCHA-1);
% %         elseif ( aZ == 1 ) % pre zeros
% %             Noise(end-N+1:end, sLCNoise(r,2), 1:numOfCHA) = Data(1:N, i:i+numOfCHA-1);
% %         elseif ( aZ == 2 ) % post zeros
% %             Noise(1:N, sLCNoise(r,2), 1:numOfCHA) = Data(1:N, i:i+numOfCHA-1);
% %         end
% %         continue;
% %     end
%     
%     % fill phase correction lines
%     r = find(s==PhaseCorrLines);
%     if ( ~isempty(r)  )
%         N = asc(i).ushSamplesInScan;
%         aZ = addPrePostZeros(asc(i));
%         if ( sLCPhaseCorr(r,2) <= maxKSpaceLineNo )
%             if ( aZ == 0 )
%                 phsCorr(1:N, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
%                     sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 1 ) % pre zeros
%                 phsCorr(end-N+1:end, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
%                     sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 2 ) % post zeros
%                 phsCorr(1:N, sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
%                     sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             end
%             phsCorrASCIndex(sLCPhaseCorr(r,2), 1:numOfCHA, sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
%                     sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = i:i+numOfCHA-1;
%         end
%         
%         % check whether it is a reflected line
%         r2 = find(s==ReflectLines);
%         if ( ~isempty(r2)  )
%              reflectPhsCorr(1, sLCPhaseCorr(r,2), 1:numOfCHA, sLCReflect(r2,4), sLCReflect(r2,5), ...
%                     sLCReflect(r2,6), sLCReflect(r2,7), sLCReflect(r2,8), sLCReflect(r2,9), sLCReflect(r2,10), sLCReflect(r2,11)) = 1;
%         end
%         
%         continue;
%     end
%     
%     % fill other lines
%     r = find(s==OtherLines);
%     if ( ~isempty(r)  )
%         N = asc(i).ushSamplesInScan;
%         aZ = addPrePostZeros(asc(i));
%         r = r(1);
%         if ( sLCOther(r,2) <= maxKSpaceLineNo )
%                        
%             if ( aZ == 0 )
%                 other(1:N, sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
%                     sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 1 ) % pre zeros
%                 other(end-N+1:end, sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
%                     sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 2 ) % post zeros
%                 other(1:N, sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
%                     sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             end
%             
%             otherASCIndex(sLCOther(r,2), 1:numOfCHA, sLCOther(r,4), sLCOther(r,5), ...
%                     sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = i:i+numOfCHA-1;
%         end
%         
%         continue;
%     end
%     
%     % fill seperate reference lines
%     r = find(s==RefLines);
%     if ( ~isempty(r)  )
%         N = asc(i).ushSamplesInScan;
%         aZ = addPrePostZeros(asc(i));
%         r = r(1);
%         if ( sLCRef(r,2) <= maxKSpaceLineNo )
%                         
%             if ( aZ == 0 )
%                 ref(1:N, sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 1 ) % pre zeros
%                 ref(end-N+1:end, sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 2 ) % post zeros
%                 ref(1:N, sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i:i+numOfCHA-1);
%             end
%             
%             refASCIndex(sLCRef(r,2), 1:numOfCHA, sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = i:i+numOfCHA-1;
%         end
%         
%         % check whether it is a reflected line
%         r2 = find(s==ReflectRefLines);
%         if ( ~isempty(r2)  )
%              reflectRef(1, sLCReflectRef(r2,2), 1:numOfCHA, sLCReflectRef(r2,4), sLCReflectRef(r2,5), ...
%                     sLCReflectRef(r2,6), sLCReflectRef(r2,7), sLCReflectRef(r2,8), sLCReflectRef(r2,9), sLCReflectRef(r2,10), sLCReflectRef(r2,11)) = 1;
%         end
%         
%         pr = find(s==RefAndImgLines);
%         if ( isempty(pr)  )
%             continue;
%         end
%     end 
%     
%     % fill the data lines
%     l = find(s==DataLines);
%     if ( ~isempty(l) )
%         N = asc(i).ushSamplesInScan;
%         
%         if ( SamplesInScan < N ) SamplesInScan = N; end
%         if ( KSpaceCentreColumn < asc(i).ushKSpaceCentreColumn ) KSpaceCentreColumn = asc(i).ushKSpaceCentreColumn; end
%         if ( KSpaceCentreLineNo < asc(i).ushKSpaceCentreLineNo ) KSpaceCentreLineNo = asc(i).ushKSpaceCentreLineNo; end
%         if ( KSpaceCentrePartitionNo < asc(i).ushKSpaceCentrePartitionNo ) KSpaceCentrePartitionNo = asc(i).ushKSpaceCentrePartitionNo; end
%         
%         aZ = addPrePostZeros(asc(i));
%         if ( sLC(l,2) <= maxKSpaceLineNo )
%             if ( aZ == 0 )
%                 kspace(1:N, sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
%                     sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 1 ) % pre zeros
%                 kspace(end-N+1:end, sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
%                     sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i:i+numOfCHA-1);
%             elseif ( aZ == 2 ) % post zeros
%                 kspace(1:N, sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
%                     sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i:i+numOfCHA-1);
%             end
%             
%             kspaceASCIndex(sLC(l,2), 1:numOfCHA, sLC(l,4), sLC(l,5), ...
%                     sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = i:i+numOfCHA-1;
%         end
%         
%         % check whether it is a reflected line
%         r2 = find(s==ReflectLines);
%         if ( ~isempty(r2)  )
%              reflect(1, sLCReflect(r2,2), 1:numOfCHA, sLCReflect(r2,4), sLCReflect(r2,5), ...
%                     sLCReflect(r2,6), sLCReflect(r2,7), sLCReflect(r2,8), sLCReflect(r2,9), sLCReflect(r2,10), sLCReflect(r2,11)) = 1;
%         end
%     end
% end

% if ( ~isempty(ref) )
%     ref = reshape(ref, sRef);
% end

toc
    % --------------------------------------------
    % parse the Eval Info Mask
    function bValue = IsNoiseScan(evalInfoMask)
        % const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4
        p = dec2bin(evalInfoMask(1));
        len = numel(p);
        if ( len >= 26 )
            if ( p(end-25) == '1' )
                bValue = 1;
            else
                bValue = 0;
            end
        else
            bValue = 0;
        end
    end

    % --------------------------------------------
    % the seperate reference lines
    function bValue = IsSeperateRef(evalInfoMask)
        % const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
        p = dec2bin(evalInfoMask(1));
        len = numel(p);
        if ( len >= 23 )
            if ( p(end-22) == '1' )
                bValue = 1;
            else
                bValue = 0;
            end
        else
            bValue = 0;
        end
    end
    % --------------------------------------------
    % parse the Eval Info Mask
    function [bNoiseLine, bSeperateRef, bRefAndImg, bPhaseCorr, bReflect, bLastScanInSlice, bRTFeedBack] = parseEvalInfoMask(evalInfoMask)
        % const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4
        % const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
        % const MdhBitField MDH_PATREFANDIMASCAN  (23);      // additonal scan for PAT reference line/partition that is also used as image scan
        % const MdhBitField MDH_REFLECT           (24);      // reflect line
        % const MdhBitField MDH_PHASCOR           (21);      // phase correction data
        % const MdhBitField MDH_LASTSCANINSLICE   (29UL);      ///< indicates  last scan in slice (needed for time stamps)
        infoMask = ['00000000000000000000000000000000'];
        p = dec2bin(evalInfoMask(1));
        len = numel(p);
        
        if ( len < 32 )
            infoMask(end-len+1:end) = p;
        else
            infoMask = p;
        end

        bNoiseLine = 0;
        if ( infoMask(end-25) == '1' ) % noise line
            bNoiseLine = 1;
        end
        
        bRTFeedBack = 0;
        if ( infoMask(end-1) == '1' ) % real time feedback, navigator probably
            bRTFeedBack = 1;
        end
               
        bSeperateRef = 0;
        if ( infoMask(end-22)=='1' ) % reference line
            bSeperateRef = 1;
        end
        
        bRefAndImg = 0;
        if ( infoMask(end-23)=='1' ) % reference and data line
            bRefAndImg = 1;
        end
        
        bPhaseCorr = 0;
        if ( infoMask(end-21)=='1' ) % phase correction line
            bPhaseCorr = 1;
        end
       
        bReflect = 0;
        if ( infoMask(end-24)=='1' ) % reflected line
            bReflect = 1;
        end
        
        bLastScanInSlice = 0;
        if ( infoMask(end-29)=='1' ) % last scan in slice
            bLastScanInSlice = 1;
        end
    end
    % --------------------------------------------
    % the acq end line
    function bValue = IsAcqEnd(evalInfoMask)
        % const MdhBitField MDH_ACQEND            ((unsigned long)0);
        bValue = (evalInfoMask(1)==0);
    end
    % --------------------------------------------
    % allocate kspace
    function [kspaceAllocated, kspaceLineASCIndex] = allocateKSpace(sLC, maxKSpaceLineNo)
        maxCol = max(sLC(:,1)); % starting from 1
        % maxLine = min(max(sLC(:,2)), maxKSpaceLineNo);
        maxLine = max(max(sLC(:,2)), maxKSpaceLineNo);
        maxCha = max(sLC(:,3));
        maxAcq = max(sLC(:,4));
        maxSlice = max(sLC(:,5));
        maxPar = max(sLC(:,6));
        maxEcho = max(sLC(:,7));
        maxPhs = max(sLC(:,8));
        maxRep = max(sLC(:,9));
        maxSet = max(sLC(:,10));
        maxSeg = max(sLC(:,11));
        % maxSeg = 1;

        % allocate the data 
        dim = [maxCol, maxLine, maxCha, maxAcq, maxSlice, maxPar, maxEcho, maxPhs, maxRep, maxSet, maxSeg]
        % kspaceAllocated = complex(single(zeros(dim)), single(zeros(dim)));
        kspaceAllocated = zeros(maxCol, maxLine, maxCha, maxAcq, maxSlice, maxPar, maxEcho, maxPhs, maxRep, maxSet, maxSeg, 'single');
        kspaceLineASCIndex = zeros(maxLine, maxCha, maxAcq, maxSlice, maxPar, maxEcho, maxPhs, maxRep, maxSet, maxSeg, 'single');
    end    
    % --------------------------------------------
    % fill the sLC
    function aSLC = fillSLC(aAsc)
        
        if ( aAsc.ushKSpaceCentreColumn == 0 )
            aSLC(1) = aAsc.ushSamplesInScan;
        else
            if ( 2*aAsc.ushKSpaceCentreColumn >= aAsc.ushSamplesInScan )
                aSLC(1) = 2*aAsc.ushKSpaceCentreColumn;
            else
                aSLC(1) = 2*(aAsc.ushSamplesInScan - aAsc.ushKSpaceCentreColumn);
            end
        end
%         aSLC(1) = max([2*aAsc.ushKSpaceCentreColumn aAsc.ushSamplesInScan]);
        aSLC(2) = aAsc.sLC(1)+1; % Line
        aSLC(3) = aAsc.ushChannelId+1; % channel
        aSLC(4:11) = aAsc.sLC(2:9)+1; % other dimensions
    end
    % --------------------------------------------
    % add pre or post zeros
    function aZ = addPrePostZeros(aAsc)
        % aZ = 1 : pre zeros
        % aZ = 2 : post zeros
        % aZ = 0 : no zeros
        if ( 2*aAsc.ushKSpaceCentreColumn == aAsc.ushSamplesInScan )
            aZ = 0;
            return;
        end
        
        if ( 2*aAsc.ushKSpaceCentreColumn < aAsc.ushSamplesInScan )
            aZ = 1;
            return;
        end
        
        if ( 2*aAsc.ushKSpaceCentreColumn > aAsc.ushSamplesInScan )
            aZ = 2;
            return;
        end
        
    end

end

