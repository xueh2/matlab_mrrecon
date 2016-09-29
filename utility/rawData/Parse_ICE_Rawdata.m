
function [kspace, kspaceASCIndex, Noise, ref, refASCIndex, phsCorr, reflect, reflectPhsCorr] = Parse_ICE_Rawdata(Data, asc) 
% read the Ice meas data
% /*--------------------------------------------------------------------------*/
% /* Definition of loop counter structure                                     */
% /* Note: any changes of this structure affect the corresponding swapping    */
% /*       method of the measurement data header proxy class (MdhProxy)       */
% /*--------------------------------------------------------------------------*/
% for the VB line
%
% typedef struct
% {
%   PACKED_MEMBER( uint16_t,  ushLine         ); /* line index                   */
%   PACKED_MEMBER( uint16_t,  ushAcquisition  ); /* acquisition index            */
%   PACKED_MEMBER( uint16_t,  ushSlice        ); /* slice index                  */
%   PACKED_MEMBER( uint16_t,  ushPartition    ); /* partition index              */
%   PACKED_MEMBER( uint16_t,  ushEcho         ); /* echo index                   */
%   PACKED_MEMBER( uint16_t,  ushPhase        ); /* phase index                  */
%   PACKED_MEMBER( uint16_t,  ushRepetition   ); /* measurement repeat index     */
%   PACKED_MEMBER( uint16_t,  ushSet          ); /* set index                    */
%   PACKED_MEMBER( uint16_t,  ushSeg          ); /* segment index  (for TSE)     */
%   PACKED_MEMBER( uint16_t,  ushIda          ); /* IceDimension a index         */
%   PACKED_MEMBER( uint16_t,  ushIdb          ); /* IceDimension b index         */
%   PACKED_MEMBER( uint16_t,  ushIdc          ); /* IceDimension c index         */
%   PACKED_MEMBER( uint16_t,  ushIdd          ); /* IceDimension d index         */
%   PACKED_MEMBER( uint16_t,  ushIde          ); /* IceDimension e index         */
% } sLoopCounter;
% /*--------------------------------------------------------------------------*/
% /*  Definition of measurement data header                                   */
% /*--------------------------------------------------------------------------*/
% typedef struct
% {
%   PACKED_MEMBER( uint32_t,     ulFlagsAndDMALength           );    // bit  0..24: DMA length [bytes]
%                                                                    // bit     25: pack bit
%                                                                    // bit 26..31: pci_rx enable flags                   4 byte
%   PACKED_MEMBER( int32_t,      lMeasUID                      );    // measurement user ID                               4
%   PACKED_MEMBER( uint32_t,     ulScanCounter                 );    // scan counter [1...]                               4
%   PACKED_MEMBER( uint32_t,     ulTimeStamp                   );    // time stamp [2.5 ms ticks since 00:00]             4
%   PACKED_MEMBER( uint32_t,     ulPMUTimeStamp                );    // PMU time stamp [2.5 ms ticks since last trigger]  4
%   PACKED_MEMBER( uint32_t,     aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK] ); // evaluation info mask field           8
%   PACKED_MEMBER( uint16_t,     ushSamplesInScan              );    // # of samples acquired in scan                     2
%   PACKED_MEMBER( uint16_t,     ushUsedChannels               );    // # of channels used in scan                        2   =32
%   sLoopCounter sLC;																					 			 // loop counters                                    28   =60
%   sCutOffData sCutOff;    																				 // cut-off values                                    4
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentreColumn         );    // centre of echo                                    2
%   PACKED_MEMBER( uint16_t,     ushCoilSelect                 );    // Bit 0..3: CoilSelect                              2
%   PACKED_MEMBER( float,        fReadOutOffcentre             );    // ReadOut offcenter value                           4
%   PACKED_MEMBER( uint32_t,     ulTimeSinceLastRF             );    // Sequence time stamp since last RF pulse           4
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentreLineNo         );    // number of K-space centre line                     2
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentrePartitionNo    );    // number of K-space centre partition                2
%   PACKED_MEMBER( uint16_t,     aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA] ); // free parameter for IceProgram   8   =88
%   PACKED_MEMBER( uint16_t,     aushFreePara[MDH_FREEHDRPARA] );    // free parameter                          4 * 2 =   8
%   sSliceData sSD;																									 // Slice Data                                       28   =124
%   PACKED_MEMBER( uint16_t,	   ushChannelId                  );    // channel Id must be the last parameter             2
%   PACKED_MEMBER( uint16_t,	   ushPTABPosNeg                 );    // negative, absolute PTAB position in [0.1 mm]      2
%                                                                    // (automatically set by PCI_TX firmware)
% } sMDH;                                                            // total length: 32 * 32 Bit (128 Byte)            128
%
% sLC: ICE loopcounter: [Line Acq Slice Partition Echo Phase Rep Set Seg Ida Idb Idc Idd Ide]
%
% EvalInfoMask
% /*--------------------------------------------------------------------------*/
% /*  Definition of EvalInfoMask:                                             */
% /*--------------------------------------------------------------------------*/
% const MdhBitField MDH_ACQEND            ((unsigned long)0);
% const MdhBitField MDH_RTFEEDBACK        (1);
% const MdhBitField MDH_HPFEEDBACK        (2);
% const MdhBitField MDH_ONLINE            (3);
% const MdhBitField MDH_OFFLINE           (4);
% const MdhBitField MDH_SYNCDATA          (5);       // readout contains synchroneous data
% const MdhBitField MDH_LASTSCANINCONCAT  (8);       // Flag for last scan in concatination
% 
% const MdhBitField MDH_RAWDATACORRECTION (10);      // Correct the rawadata with the rawdata correction factor
% const MdhBitField MDH_LASTSCANINMEAS    (11);      // Flag for last scan in measurement
% const MdhBitField MDH_SCANSCALEFACTOR   (12);      // Flag for scan specific additional scale factor
% const MdhBitField MDH_2NDHADAMARPULSE   (13);      // 2nd RF exitation of HADAMAR
% const MdhBitField MDH_REFPHASESTABSCAN  (14);      // reference phase stabilization scan         
% const MdhBitField MDH_PHASESTABSCAN     (15);      // phase stabilization scan
% const MdhBitField MDH_D3FFT             (16);      // execute 3D FFT         
% const MdhBitField MDH_SIGNREV           (17);      // sign reversal
% const MdhBitField MDH_PHASEFFT          (18);      // execute phase fft     
% const MdhBitField MDH_SWAPPED           (19);      // swapped phase/readout direction
% const MdhBitField MDH_POSTSHAREDLINE    (20);      // shared line               
% const MdhBitField MDH_PHASCOR           (21);      // phase correction data    
% const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
% const MdhBitField MDH_PATREFANDIMASCAN  (23);      // additonal scan for PAT reference line/partition that is also used as image scan
% const MdhBitField MDH_REFLECT           (24);      // reflect line              
% const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4        
% const MdhBitField MDH_SHARENOW          (26);      // all lines are acquired from the actual and previous e.g. phases
% const MdhBitField MDH_LASTMEASUREDLINE  (27);      // indicates that the current line is the last measured line of all succeeding e.g. phases
% const MdhBitField MDH_FIRSTSCANINSLICE  (28);      // indicates first scan in slice (needed for time stamps)
% const MdhBitField MDH_LASTSCANINSLICE   (29);      // indicates  last scan in slice (needed for time stamps)
% const MdhBitField MDH_TREFFECTIVEBEGIN  (30);      // indicates the begin time stamp for TReff (triggered measurement)
% const MdhBitField MDH_TREFFECTIVEEND    (31);      // indicates the   end time stamp for TReff (triggered measurement)
% const MdhBitField MDH_MDS_REF_POSITION  (32);      // indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
% const MdhBitField MDH_SLC_AVERAGED      (33);      // indicates avveraged slice for slice partial averaging scheme
% 
% const MdhBitField MDH_FIRST_SCAN_IN_BLADE       (40);  // Marks the first line of a blade
% const MdhBitField MDH_LAST_SCAN_IN_BLADE        (41);  // Marks the last line of a blade
% const MdhBitField MDH_LAST_BLADE_IN_TR          (42);  // Set for all lines of the last BLADE in each TR interval
% 
% const MdhBitField MDH_RETRO_LASTPHASE           (45);  // Marks the last phase in a heartbeat
% const MdhBitField MDH_RETRO_ENDOFMEAS           (46);  // Marks an ADC at the end of the measurement
% const MdhBitField MDH_RETRO_REPEATTHISHEARTBEAT (47);  // Repeat the current heartbeat when this bit is found
% const MdhBitField MDH_RETRO_REPEATPREVHEARTBEAT (48);  // Repeat the previous heartbeat when this bit is found
% const MdhBitField MDH_RETRO_ABORTSCANNOW        (49);  // Just abort everything
% const MdhBitField MDH_RETRO_LASTHEARTBEAT       (50);  // This adc is from the last heartbeat (a dummy)
% const MdhBitField MDH_RETRO_DUMMYSCAN           (51);  // This adc is just a dummy scan, throw it away
% const MdhBitField MDH_RETRO_ARRDETDISABLED      (52);  // Disable all arrhythmia detection when this bit is found
% 
% //-----------------------------------------------------------------------------
% // Definition of EvalInfoMask for COP:
% //-----------------------------------------------------------------------------
% #define	MDH_COP_ACQEND            (0x00000001)
% #define MDH_COP_RTFEEDBACK        (0x00000002)
% #define MDH_COP_HPFEEDBACK        (0x00000004)
% #define	MDH_COP_ONLINE            (0x00000008)
% #define	MDH_COP_OFFLINE           (0x00000010)
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


tic
S0 = size(Data);
numOfLines = S0(2);

% find the noise lines and seperate ref lines if any
numOfNoiseLines = 0; 
NoiseLines = zeros(numOfLines,1, 'int32');
sLCNoise = zeros(numOfLines, 11, 'int32');

numOfSeperateRefLines = 0; 
RefLines = zeros(numOfLines,1, 'int32');
sLCRef = zeros(numOfLines, 11, 'int32');

numOfPhaseCorrLines = 0; 
PhaseCorrLines = zeros(numOfLines,1, 'int32');
sLCPhaseCorr = zeros(numOfLines, 11, 'int32');

numOfReflectLines = 0; 
ReflectLines = zeros(numOfLines,1, 'int32');
sLCReflect = zeros(numOfLines, 11, 'int32');

numOfDataLines = 0; 
DataLines = zeros(numOfLines,1, 'int32');
sLC = zeros(numOfLines, 11, 'int32');

kspaceCentreLineNo = zeros(numOfLines,1, 'int32');

for i=1:numOfLines
      
    if ( IsAcqEnd(asc(i).aulEvalInfoMask(1)) )
        continue;
    end
    
    kspaceCentreLineNo(i) = asc(i).ushKSpaceCentreLineNo;
    
    [bNoiseLine, bSeperateRef, bPhaseCorr, bReflect] = parseEvalInfoMask(asc(i).aulEvalInfoMask(1));
    
    if ( bNoiseLine  )
        numOfNoiseLines = numOfNoiseLines + 1;
        NoiseLines(numOfNoiseLines) = i;
        sLCNoise(numOfNoiseLines,:) = fillSLC(asc(i));
        continue;
    end
    
    if ( bSeperateRef )
        numOfSeperateRefLines = numOfSeperateRefLines + 1;
        RefLines(numOfSeperateRefLines) = i;
        sLCRef(numOfSeperateRefLines,:) = fillSLC(asc(i));
        continue;
    end 

    if ( bPhaseCorr )
        numOfPhaseCorrLines = numOfPhaseCorrLines + 1;
        PhaseCorrLines(numOfPhaseCorrLines) = i;
        sLCPhaseCorr(numOfPhaseCorrLines,:) = fillSLC(asc(i));
        if ( bReflect ) % if it is a reflected phase correction lines
            numOfReflectLines = numOfReflectLines + 1;
            ReflectLines(numOfReflectLines) = i;
            sLCReflect(numOfReflectLines,:) = fillSLC(asc(i));
        end
        continue;
    end 
    
    if ( bReflect ) % if it is a reflected data lines
        numOfReflectLines = numOfReflectLines + 1;
        ReflectLines(numOfReflectLines) = i;
        sLCReflect(numOfReflectLines,:) = fillSLC(asc(i));
    end 
    
    numOfDataLines = numOfDataLines + 1;
    DataLines(numOfDataLines) = i;
    sLC(numOfDataLines,:) = fillSLC(asc(i));
end

sLCNoise = sLCNoise(1:numOfNoiseLines, :); NoiseLines = NoiseLines(1:numOfNoiseLines);

sLCRef = sLCRef(1:numOfSeperateRefLines, :); RefLines = RefLines(1:numOfSeperateRefLines);

sLCPhaseCorr = sLCPhaseCorr(1:numOfPhaseCorrLines, :); PhaseCorrLines = PhaseCorrLines(1:numOfPhaseCorrLines);

sLCReflect = sLCReflect(1:numOfReflectLines, :); ReflectLines = ReflectLines(1:numOfReflectLines);

sLC = sLC(1:numOfDataLines, :); DataLines = DataLines(1:numOfDataLines);

maxKSpaceLineNo = 2*max(kspaceCentreLineNo);
[kspace, kspaceASCIndex] = allocateKSpace(sLC, maxKSpaceLineNo);

% SEG = size(kspace, 11);

% [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
% 11 dimension array

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
    [ref, refASCIndex] = allocateKSpace(sLCRef, maxKSpaceLineNo);
%     sRef = size(ref);
%     NRef = prod(sRef(2:11));
%     ref = reshape(ref, [sRef(1) NRef]);    
else
    ref = [];
    refASCIndex = [];
end

% phase correction lines
if ( numOfPhaseCorrLines > 0 )
    [phsCorr, phsCorrASCIndex] = allocateKSpace(sLCPhaseCorr, maxKSpaceLineNo);
    ss = size(phsCorr);
    reflectPhsCorr = zeros([1 ss(2:end)], 'int8');
else
    phsCorr = [];
    phsCorrASCIndex = [];
    reflectPhsCorr = [];
end

% phase correction lines
if ( numOfReflectLines > 0 )
    ss = size(kspace);
    reflect = zeros([1 ss(2:end)], 'int8');
else
    reflect = [];
end

% put the data into the array

l = 0;
r = 0;
n = 0;
for i=1:numOfLines
    
    if ( mod(i, 1e3) == 0 )
        disp([num2str(i) ' lines parsed ... ']);
    end
    
    if ( IsAcqEnd(asc(i).aulEvalInfoMask(1)) )
        continue;
    end

    % fill noise lines
    r = find(i==NoiseLines);
    if ( ~isempty(r) )
        n = n + 1;
        N = asc(i).ushSamplesInScan;
        aZ = addPrePostZeros(asc(i));
        if ( aZ == 0 )
            Noise(1:N, sLCNoise(r,2), sLCNoise(r,3)) = Data(1:N, i);
        elseif ( aZ == 1 ) % pre zeros
            Noise(end-N+1:end, sLCNoise(r,2), sLCNoise(r,3)) = Data(1:N, i);
        elseif ( aZ == 2 ) % post zeros
            Noise(1:N, sLCNoise(r,2), sLCNoise(r,3)) = Data(1:N, i);
        end
        continue;
    end
    
    % fill seperate reference lines
    r = find(i==RefLines);
    if ( ~isempty(r)  )
        i;
        N = asc(i).ushSamplesInScan;
        aZ = addPrePostZeros(asc(i));
        r = r(1);
        if ( sLCRef(r,2) <= maxKSpaceLineNo )
            
%             indRef = sub2ind(sRef(2:11), sLCRef(r,2), sLCRef(r,3), sLCRef(r,4), sLCRef(r,5), ...
%                     sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11));
            
            if ( aZ == 0 )
                ref(1:N, sLCRef(r,2), sLCRef(r,3), sLCRef(r,4), sLCRef(r,5), ...
                    sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i);
%                 ref(1:N, indRef) = Data(1:N, i);
            elseif ( aZ == 1 ) % pre zeros
                ref(end-N+1:end, sLCRef(r,2), sLCRef(r,3), sLCRef(r,4), sLCRef(r,5), ...
                    sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i);
            elseif ( aZ == 2 ) % post zeros
                ref(1:N, sLCRef(r,2), sLCRef(r,3), sLCRef(r,4), sLCRef(r,5), ...
                    sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i);
            end
            
            refASCIndex(sLCRef(r,2), sLCRef(r,3), sLCRef(r,4), sLCRef(r,5), ...
                    sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = i;
        end
        
        continue;
    end
    
    % fill phase correction lines
    r = find(i==PhaseCorrLines);
    if ( ~isempty(r)  )
        N = asc(i).ushSamplesInScan;
        aZ = addPrePostZeros(asc(i));
        if ( sLCPhaseCorr(r,2) <= maxKSpaceLineNo )
            if ( aZ == 0 )
                phsCorr(1:N, sLCPhaseCorr(r,2), sLCPhaseCorr(r,3), sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i);
            elseif ( aZ == 1 ) % pre zeros
                phsCorr(end-N+1:end, sLCPhaseCorr(r,2), sLCPhaseCorr(r,3), sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i);
            elseif ( aZ == 2 ) % post zeros
                phsCorr(1:N, sLCPhaseCorr(r,2), sLCPhaseCorr(r,3), sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = Data(1:N, i);
            end
            phsCorrASCIndex(sLCPhaseCorr(r,2), sLCPhaseCorr(r,3), sLCPhaseCorr(r,4), sLCPhaseCorr(r,5), ...
                    sLCPhaseCorr(r,6), sLCPhaseCorr(r,7), sLCPhaseCorr(r,8), sLCPhaseCorr(r,9), sLCPhaseCorr(r,10), sLCPhaseCorr(r,11)) = i;
        end
        
        % check whether it is a reflected line
        r2 = find(i==ReflectLines);
        if ( ~isempty(r2)  )
             reflectPhsCorr(1, sLCReflect(r2,2), sLCReflect(r2,3), sLCReflect(r2,4), sLCReflect(r2,5), ...
                    sLCReflect(r2,6), sLCReflect(r2,7), sLCReflect(r2,8), sLCReflect(r2,9), sLCReflect(r2,10), sLCReflect(r2,11)) = 1;
        end
        
        continue;
    end
    
    % fill the data lines
    l = find(i==DataLines);
    if ( ~isempty(l) )
        N = asc(i).ushSamplesInScan;
        aZ = addPrePostZeros(asc(i));
        if ( sLC(l,2) <= maxKSpaceLineNo )
            if ( aZ == 0 )
                kspace(1:N, sLC(l,2), sLC(l,3), sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i);
            elseif ( aZ == 1 ) % pre zeros
                kspace(end-N+1:end, sLC(l,2), sLC(l,3), sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i);
            elseif ( aZ == 2 ) % post zeros
                kspace(1:N, sLC(l,2), sLC(l,3), sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = Data(1:N, i);
            end
            
            kspaceASCIndex(sLC(l,2), sLC(l,3), sLC(l,4), sLC(l,5), ...
                    sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = i;
        end
        
        % check whether it is a reflected line
        r2 = find(i==ReflectLines);
        if ( ~isempty(r2)  )
             reflect(1, sLCReflect(r2,2), sLCReflect(r2,3), sLCReflect(r2,4), sLCReflect(r2,5), ...
                    sLCReflect(r2,6), sLCReflect(r2,7), sLCReflect(r2,8), sLCReflect(r2,9), sLCReflect(r2,10), sLCReflect(r2,11)) = 1;
        end
    end
end

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
    function [bNoiseLine, bSeperateRef, bPhaseCorr, bReflect] = parseEvalInfoMask(evalInfoMask)
        % const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4
        % const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
        % const MdhBitField MDH_PATREFANDIMASCAN  (23);      // additonal scan for PAT reference line/partition that is also used as image scan
        % const MdhBitField MDH_REFLECT           (24);      // reflect line
        % const MdhBitField MDH_PHASCOR           (21);      // phase correction data
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
        
        bSeperateRef = 0;
        if ( infoMask(end-22)=='1' & infoMask(end-23)~='1' ) % reference line
            bSeperateRef = 1;
        end
        
        bPhaseCorr = 0;
        if ( infoMask(end-21)=='1' ) % phase correction line
            bPhaseCorr = 1;
        end
       
        bReflect = 0;
        if ( infoMask(end-24)=='1' ) % reflected line
            bReflect = 1;
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

