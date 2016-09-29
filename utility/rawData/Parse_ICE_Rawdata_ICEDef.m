
%% VB17
% /*---------------------------------------------------------------------------*/
% /*  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential        */
% /*---------------------------------------------------------------------------*/
% /*
%  * Project: NUMARIS/4
%  *    File: \n4\pkg\MrServers\MrProtSrv\MrProt\SeqDefines.h
%  * Version:
%  *  Author: Comp_ProBe
%  *    Date: n.a.
%  *
%  *    Lang: CPP
%  *
%  * Descrip:
%  *
%  * Functns: n.a.
%  *
%  *---------------------------------------------------------------------------*/
% 
% 
% /*--------------------------------------------------------------------------*/
% /* Include control                                                          */
% /*--------------------------------------------------------------------------*/
% #ifndef __SeqDefines
% #define __SeqDefines
% 
% //-----------------------------------------------------------------------------
% // Class prototypes and local typedefs for it.
% // (see also the NUCLEUS_xxx defines in this file)
% //-----------------------------------------------------------------------------
% #ifdef __cplusplus
%   class   MeasNucleus;
%   class   MeasNucleusSet;
%   class   MeasKnownNuclei;
%   typedef MeasNucleus      SEQ_MEASNUCLEUS;
%   typedef MeasNucleusSet   SEQ_MEASNUCLEUSSET;
%   typedef MeasKnownNuclei  SEQ_MEASKNOWNNUCLEI;
% #endif
% 
% 
% 
% //-----------------------------------------------------------------------------
% // Structure to provide the "namespace" SEQ
% //-----------------------------------------------------------------------------
% // stupid stabs anna cant handle namespaces (the complete tool is shitty)
% 
% struct SEQ
% {
%   //---------------------------------------------------------------------------
%   //  Enums for sequence parameter
%   //---------------------------------------------------------------------------
%   enum GlobalRFactor
%   {
%     GLOBALRFACTOR_DEFAULT = 0x01,
%     GLOBALRFACTOR_1       = 0x02,
%     GLOBALRFACTOR_2       = 0x04,
%     GLOBALRFACTOR_3       = 0x08
%   };
% 
%   enum Switch
%   {
%     NO      = 0x1,
%     OFF     = 0x1,
%     DISABLE = 0x1,
%     YES     = 0x2,
%     ON      = 0x2,
%     ENABLE  = 0x2,
%     ALLOWED = 0x2
%   };
% 
%   //--------------------------------------------------------------------------------------
%   // Allows the sequence to control usage of a parameter the value of which may
%   // be defined by the system (e.g. adjustments and/or sequence; i.e. AUTO mode ) 
%   // or may be set to a fix value by the user and/or the sequence (e.g. by init settings).
%   // The limits and values of the parameter value itself is controlled by separate
%   // SeqLim and MrProt entries.
%   // To control the AUTO mode behaviour in SeqLim and MrProt use the 
%   // xxxValid() methods for this parameter.
%   // e.g.:
%   //  if ( adjustWithToleranceValid() ) { /* value for adjustWithTolerance() was defined by user or sequence */ }
%   //  if ( !xxxValid() )                { /* value for xxx() is in AUTO mode and may have system specific values in different use cases */ }
%   //--------------------------------------------------------------------------------------
%   enum ParamValidFlag
%   {
%     PARAM_AUTO      = 0x1,
%     PARAM_VALID     = 0x2
%   };
% 
%   enum DisplayMode
%   {
%     DM_OFF  = 0x01,
%     DM_SHOW = 0x02,
%     DM_EDIT = 0x04
%   };
% 
%   enum UsedADC
%   {
%     ADC1  = 0x00000001,
%     ADC2  = 0x00000002,
%     ADC3  = 0x00000004,
%     ADC4  = 0x00000008,
%     ADC5  = 0x00000010,
%     ADC6  = 0x00000020,
%     ADC7  = 0x00000040,
%     ADC8  = 0x00000080,
%     ADC9  = 0x00000100,
%     ADC10 = 0x00000200,
%     ADC11 = 0x00000400,
%     ADC12 = 0x00000800,
%     ADC13 = 0x00001000,
%     ADC14 = 0x00002000,
%     ADC15 = 0x00004000,
%     ADC16 = 0x00008000,
%     ADC17 = 0x00010000,
%     ADC18 = 0x00020000,
%     ADC19 = 0x00040000,
%     ADC20 = 0x00080000,
%     ADC21 = 0x00100000,
%     ADC22 = 0x00200000,
%     ADC23 = 0x00400000,
%     ADC24 = 0x00800000,
%     ADC25 = 0x01008000,
%     ADC26 = 0x02000000,
%     ADC27 = 0x04000000,
%     ADC28 = 0x08000000,
%     ADC29 = 0x10000000,
%     ADC30 = 0x20000000,
%     ADC31 = 0x40000000,
%     ADC32 = 0x80000000
%   };
% 
%   enum AveragingMode
%   {
%     INNER_LOOP = 0x01,
%     OUTER_LOOP = 0x02
%   };
% 
%   enum Dimension
%   {
%     DIM_1 = 0x01,
%     DIM_2 = 0x02,
%     DIM_3 = 0x04
%   };
% 
%   enum MultiSliceMode
%   {
%     MSM_SEQUENTIAL  = 0x01,
%     MSM_INTERLEAVED = 0x02,
%     MSM_SINGLESHOT  = 0x04
%   };
% 
%   enum SimultaneousExcitation
%   {
%     HADAMARD_OFF     = 0x01,
%     HADAMARD_2SLICES = 0x02
%   };
% 
%   enum PartialFourierFactor
%   {
%     PF_HALF = 0x01,
%     PF_5_8  = 0x02,
%     PF_6_8  = 0x04,
%     PF_7_8  = 0x08,
%     PF_OFF  = 0x10,
%     PF_AUTO = 0x20
%   };
% 
%   enum POCS //(projection onto convex sets recon mode)
%   {
%     POCS_OFF        = 0x01,
%     POCS_READ_SLICE = 0x02,
%     POCS_READ_PHASE = 0x04
%   };
% 
%   enum FatSuppression
%   {
%     FAT_SATURATION          = 0x01,
%     WATER_EXCITATION        = 0x02,
%     FAT_SUPPRESSION_OFF     = 0x04,
%     FAT_SATURATION_QUICK    = 0x08,
%     WATER_EXCITATION_FAST   = 0x10,
%     FAT_SUPPRESSION_OPTIMAL = 0x20
%   };
% 
%   enum FatSatMode
%   {
%     FAT_SAT_WEAK         = 0x01,
%     FAT_SAT_STRONG       = 0x02
%   };
% 
%   enum AutoAlignAlgo
%   {
%     AA_ALGO_NONE         = 0x01,
%     AA_HEAD              = 0x02,
%     AA_KNEE              = 0x04,
%     AA_HEAD_ATLAS        = 0x08
%   };
% 
%   enum AutoAlignRegion
%   { 
%     AA_REGION_NONE       = 0x01,
%     AA_REGION_HEAD       = 0x02,
%     AA_REGION_SPINE      = 0x04,
%     AA_REGION_KNEE       = 0x08
%   };
% 
%   enum AutoAlignRefMatrix
%   {
%     AA_REF_NONE          = 0x01,
%     AA_REF_SCULL         = 0x02,
%     AA_REF_BRAIN         = 0x04,
%     AA_REF_IAC           = 0x08,
%     AA_REF_ORBITA_TRA    = 0x10,
%     AA_REF_ORBITA_SAGL   = 0x20,
%     AA_REF_ORBITA_SAGR   = 0x40,
%     AA_REF_BRAINATLAS    = 0x80,
%     AA_REF_VERTEBRADISKS = 0x100,
%     AA_REF_SAG_COR_TRA   = 0x200
%   };
% 
%   enum AutoAlignMatrixTag
%   {
%     AA_MAT_NONE            = 0x01,
%     AA_MAT_H_SCULL         = 0x02,
%     AA_MAT_H_BRAIN         = 0x04,
%     AA_MAT_H_IAC           = 0x08,
%     AA_MAT_H_ORBITA_TRA    = 0x10,
%     AA_MAT_H_ORBITA_SAGL   = 0x20,
%     AA_MAT_H_ORBITA_SAGR   = 0x40,
%     AA_MAT_H_BRAINATLAS    = 0x80,
%     AA_MAT_S_VERTEBRADISKS = 0x100,
%     AA_MAT_K_SAG_COR_TRA   = 0x200
%   };
% 
%   enum AutoAlignLicState
%   {
%     AA_LIC_AAMR            = 0x01,
%     AA_LIC_HEAD            = 0x02,
%     AA_LIC_KNEE            = 0x04,
%     AA_LIC_SPINE           = 0x08
%   };
% 
%   enum WaterSuppression
%   {
%     WATER_SATURATION          = 0x01,
%     FAT_EXCITATION            = 0x02,
%     WATER_SUPPRESSION_OFF     = 0x04,
%     WATER_SATURATION_QUICK    = 0x08,
%     WATER_SUPPRESSION_PARTIAL = 0x10,
%     WATER_SUPPRESSION_WEAK    = 0x20,
%     WATER_SUPPRESSION_RF_OFF  = 0x40
%   };
% 
%   enum Dixon
%   {
%     DIXON_NONE             = 0x01,
%     DIXON_WATER_IMAGE      = 0x02,
%     DIXON_FAT_IMAGE        = 0x04,
%     DIXON_WATER_FAT_IMAGES = 0x08
%   };
% 
%   enum Inversion
%   {
%     SLICE_SELECTIVE  = 0x01,
%     VOLUME_SELECTIVE = 0x02,
%     INVERSION_OFF    = 0x04,
%     T2_SELECTIVE		 = 0x08,
%     VOLUME_SELECTIVE_DIR = 0x10
%   };
% 
%   enum T2Prep
%   {
%     T2_PREPARATION_OFF  = 0x01,
%     T2_PREPARATION_ON   = 0x02,
%   };
% 
%   enum TIScout
%   {
%     TI_SCOUT_OFF  = 0x01,
%     TI_SCOUT_ON   = 0x02,
%   };
% 
%   enum SeriesMode
%   {
%     ASCENDING   = 0x01,
%     DESCENDING  = 0x02,
%     INTERLEAVED = 0x04, 
%     AUTOMATIC   = 0x08,
%     APEXTOBASE  = 0x10,
%     BASETOAPEX  = 0x20
%   };
% 
%   enum SliceOrientation
%   {
%     ORTHOGONAL     = 0x01,  // includes SAG/COR/TRA
%     SINGLE_OBLIQUE = 0x02,  // includes ORTHOGONAL
%     DOUBLE_OBLIQUE = 0x04,  // includes ORTHOGONAL and SINGLE_OBLIQUE
%     SAG            = 0x08,
%     COR            = 0x10,
%     TRA            = 0x20
%   };
% 
%   enum MdsSliceOrientation
%   {
%     MDS_SAG            = 0x01,
%     MDS_COR            = 0x02,
%     MDS_TRA            = 0x04
%   };
% 
%   enum RFPulseType
%   {
%     RF_FAST      = 0x01,
%     RF_NORMAL    = 0x02,
%     RF_LOW_SAR   = 0x04,
%     RF_OPTIMIZED = 0x08
%   };
% 
%   enum RFPulseShape
%   {
%     RF_BINOMIAL     = 0x01,
%     RF_HYPERSECANT  = 0x02,
%     RF_RECTANGULAR  = 0x04,
%     RF_HYPERTANGENT = 0x08
%   };
% 
%   enum Gradients
%   {
%     GRAD_FAST                  = 0x01,
%     GRAD_NORMAL                = 0x02,
%     GRAD_WHISPER               = 0x04,
%     GRAD_FAST_GSWD_RISETIME    = (GRAD_FAST    | 0x10),
%     GRAD_NORMAL_GSWD_RISETIME  = (GRAD_NORMAL  | 0x10),
%     GRAD_WHISPER_GSWD_RISETIME = (GRAD_WHISPER | 0x10)
%   };
% 
%   enum GradientAxis
%   {
%     AXIS_UNDEFINED = 0x00,
%     AXIS_PHASE     = 0x01,
%     AXIS_READOUT   = 0x02,
%     AXIS_SLICE     = 0x04
%   };
% 
%   enum Direction
%   {
%     DIR_ASCENDING  = 0x01,
%     DIR_DESCENDING = 0x02
%   };
% 
%   enum ReconstructionMode
%   {
%     RECONMODE_MAGNITUDE  = 0x01,
%     RECONMODE_PHASE      = 0x02,
%     RECONMODE_REAL_PART  = 0x04,
%     RECONMODE_MAGN_PHASE = 0x08,
%     RECONMODE_REAL_PHASE = 0x10,
%     RECONMODE_PSIR       = 0x20
%   };
% 
%   enum RegriddingMode
%   {
%     REGRID_NONE        = 0x01,
%     REGRID_TRAPEZOIDAL = 0x02,
%     REGRID_SINUSOIDAL  = 0x04
%   };
% 
%   enum FilterMode
%   {
%     FILTER_WEAK   = 0x01,
%     FILTER_MEDIUM = 0x02,
%     FILTER_STRONG = 0x04,
%     FILTER_FREE   = 0x08,
%     FILTER_SHARP  = 0x10,
%     FILTER_SMOOTH = 0x20
%   };
% 
%   enum EllipticalFilterMode
%   {
%     ELLIPTICAL_FILTER_INPLANE   = 0x01,
%     ELLIPTICAL_FILTER_VOLUME    = 0x02
%   };
% 
%   enum Base
%   {
%     BASE2 = 1
%   };
% 
%   enum Increment
%   {
%     INC_NORMAL           = 1,
%     INC_BASE2            = 2,
%     INC_FIX              = 3,
%     INC_TURBO_FACTOR     = 4,
%     INC_EPI_FACTOR       = 5,
%     INC_SEGMENTED        = 6,
%     INC_TGSE_FACTOR      = 7,
%     INC_GRE_SEGMENTS     = 8,
%     INC_SINGLESHOT       = 9,
%     INC_TWICE_EPI_FACTOR = 10,
%     INC_64               = 11,
%     INC_32               = 12,
%     INC_16               = 13,
%     /// UILink / KSpace taking care of Reflines and IPAT
%     INC_NORMAL_IPAT           = 14,
%     INC_BASE2_IPAT            = 15,
%     INC_FIX_IPAT              = 16,
%     INC_TURBO_FACTOR_IPAT     = 17,
%     INC_EPI_FACTOR_IPAT       = 18,
%     INC_SEGMENTED_IPAT        = 19,
%     INC_TGSE_FACTOR_IPAT      = 20,
%     INC_GRE_SEGMENTS_IPAT     = 21,
%     INC_SINGLESHOT_IPAT       = 22,
%     INC_TWICE_EPI_FACTOR_IPAT = 23,
%     INC_64_IPAT               = 24,
%     INC_32_IPAT               = 25,
%     INC_16_IPAT               = 26,
%   };
% 
%   enum SequenceCard
%   {
%     SEQUENCE_CARD_NONE         = 1,
%     SEQUENCE_CARD_IMAGING      = 2,
%     SEQUENCE_CARD_SPECTROSCOPY = 3
%   };
% 
%   enum ApplicationCard
%   {
%     APPLICATION_CARD_NONE         = 1,
%     APPLICATION_CARD_FMRI         = 2,
%     APPLICATION_CARD_DIFF         = 3,
%     APPLICATION_CARD_TOF_CE       = 4,
%     APPLICATION_CARD_PERF         = 5,
%     APPLICATION_CARD_PC_FLOW      = 6,
%     APPLICATION_CARD_TOF_PCANGIO  = 7,
%     APPLICATION_CARD_INLINE       = 8,
%     APPLICATION_CARD_EVA_CUSTOM   = 9,
%     APPLICATION_CARD_MREPORT      = 10
%   };
% 
%   enum ApplicationCardName
%   {
%     APPLICATION_CARD_NAME_NONE       = 1,
%     APPLICATION_CARD_NAME_BOLD       = 2,
%     APPLICATION_CARD_NAME_FMRI       = 3,
%     APPLICATION_CARD_NAME_PERF       = 4,
%     APPLICATION_CARD_NAME_DIFF       = 5,
%     APPLICATION_CARD_NAME_TOF        = 6,
%     APPLICATION_CARD_NAME_PC         = 7,
%     APPLICATION_CARD_NAME_FLOW       = 8,
%     APPLICATION_CARD_NAME_CE         = 9,
%     APPLICATION_CARD_NAME_MRA        = 10,
%     APPLICATION_CARD_NAME_INLINE     = 11,
%     APPLICATION_CARD_NAME_EVA_CUSTOM = 12,
%     APPLICATION_CARD_NAME_MREPORT    = 13
%   };
% 
%   enum DiffusionMode
%   {
%     DIFFMODE_NONE             = 0x01,
%     DIFFMODE_ONE_SCAN_TRACE   = 0x02,
%     DIFFMODE_ORTHOGONAL       = 0x04,
%     DIFFMODE_SLICE            = 0x08,
%     DIFFMODE_READ             = 0x10,
%     DIFFMODE_PHASE            = 0x20,
%     DIFFMODE_THREE_SCAN_TRACE = 0x40,
%     DIFFMODE_FREE             = 0x80,
%     DIFFMODE_TENSOR           = 0x100,
%     DIFFMODE_DIAGONAL         = 0x200
%   };
% 
%   enum DiffusionScheme
%   {
%       DIFFSCHEME_BIPOLAR = 0x01,
%       DIFFSCHEME_MONOPOLAR = 0x02,
%       DIFFSCHEME_MONOPOLAR_PLUS = 0x04,
%   };
% 
%   ////////////////////////////////////////////////////////////////////////////
%   //
%   // Values for the bit-field ucBOLDParadigmArray[].
%   // For each repetition (within index range [0,1023], 1st measurement corresponds to repetition index 0)
%   // there are two element-bits representing the paradigm for the corresponding repetition.
%   //
%   // The paradigm is periodic after a number of initial measurements, this periodicity
%   // (MrProt::ParadigmPeriodicity())can therefore be used
%   // to determine the paradigm also for repetition indices > 1023.
%   //
%   // MrProt has two methods to get/set the element-bits:
%   //
%   //    SEQ::ParadigmElem MrProt::ParadigmByRepetition(long lRepetition)
%   //        returns the paradigm-element-bits
%   //        corresponding to the given repetition
%   //
%   //    SEQ::ParadigmElem MrProt::ParadigmByRepetition(long lRepetition,
%   //                                                    SEQ::ParadigmElem eParadigm)
%   //        will set the paradigm-element-bits
%   //        corresponding to the given repetition
%   //
%   // MrProt has a method to determine if the actual protocol
%   // has an active BOLD measurement
%   //    bool isBOLDMeasurement()
%   //        returns true, if lParadigmPeriodicity != 0
%   //
%   ////////////////////////////////////////////////////////////////////////////
%   enum ParadigmElem
%   {
%     PARADIGM_IGNORE,            // repetition is ignored
%     PARADIGM_ACTIVATED,         // repetition is baseline
%     PARADIGM_BASELINE           // repetition is activated
%   };
% 
% 
%   enum MainOrientation
%   {
%     SAGITTAL   = 0,
%     CORONAL    = 1,
%     TRANSVERSE = 2
%   };
% 
%   enum MainDirection
%   {
%     R_TO_L = SAGITTAL,         // (R)ight to (L)eft
%     L_TO_R = 4,                // (L)eft to (R)ight
%     A_TO_P = CORONAL,          // (A)nterior to (P)osterior
%     P_TO_A = 8,                // (P)osterior to (A)nterior
%     F_TO_H = TRANSVERSE,       // (F)eet to (H)ead
%     H_TO_F = 16                // (H)ead to (F)eet
%   };
% 
%   enum RotationCode
%   {
%     SAG_TO_TRA_TO_COR = 1,
%     SAG_TO_COR_TO_TRA = 2,
%     COR_TO_SAG_TO_TRA = 3,
%     COR_TO_TRA_TO_SAG = 4,
%     TRA_TO_COR_TO_SAG = 5,
%     TRA_TO_SAG_TO_COR = 6
%   };
% 
%   enum PhaseCyclingType
%   {
%     PHASE_CYCLING_NONE             = 0x01,
%     PHASE_CYCLING_AUTO             = 0x02,
%     PHASE_CYCLING_TWOSTEP          = 0x04,
%     PHASE_CYCLING_EIGHTSTEP        = 0x08,
%     PHASE_CYCLING_EXORCYCLE        = 0x10,
%     PHASE_CYCLING_SIXTEENSTEP_EXOR = 0x20,
%     PHASE_CYCLING_SMART_AVERAGE    = 0x40
%   };
% 
%   enum PhapsMode
%   {
%     PHAPS_NONE       = 0x01,
%     PHAPS_ON         = 0x02,
%     PHAPS_SUM_ONLY   = 0x04,
%     PHAPS_DIFF_ONLY  = 0x08
%   };
% 
%   enum PhaseEncodingType
%   {
%     PHASE_ENCODING_FULL       = 0x01,
%     PHASE_ENCODING_ELLIPTICAL = 0x02,
%     PHASE_ENCODING_WEIGHTED   = 0x04
%   };
% 
%   enum RespirationPhase
%   {
%     PHASE_INSPIRATION = 0x01,
%     PHASE_EXPIRATION  = 0x02
%   };
% 
%   enum PhysioSignal
%   {
%     SIGNAL_NONE        = 0x01,
%     SIGNAL_EKG         = 0x02,
%     SIGNAL_PULSE       = 0x04,
%     SIGNAL_EXT         = 0x08,
%     SIGNAL_CARDIAC     = 0x0E,  /* the sequence usually takes this */
%     SIGNAL_RESPIRATION = 0x10,
%     SIGNAL_ALL         = 0x1E,
%     SIGNAL_EKG_AVF     = 0x20
%   };
% 
%   enum PhysioMethod
%   {
%     METHOD_NONE        = 0x01,
%     METHOD_TRIGGERING  = 0x02,
%     METHOD_GATING      = 0x04,
%     METHOD_RETROGATING = 0x08,
%     METHOD_SOPE        = 0x10,
%     METHOD_ALL         = 0x1E
%   };
%    
%   enum PhysioNativeMode
%   {
%     NATIVE_NONE        = 0x01,
%     NATIVE_NORMAL      = 0x02,
%     NATIVE_DYNAMIC     = 0x04,
%     NATIVE_TT_SCOUT    = 0x08,
%     NATIVE_TT_AUTO     = 0x10,
%     NATIVE_3D_AUTO     = 0x20
%   };
% 
%   enum NativeFlowSensitivity
%   {
%     NATIVE_FLOW_SENSITIVITY_OFF       = 0x01,
%     NATIVE_FLOW_SENSITIVITY_WEAK      = 0x02,
%     NATIVE_FLOW_SENSITIVITY_MEDIUM    = 0x04,
%     NATIVE_FLOW_SENSITIVITY_STRONG    = 0x08        
%   };
% 
%   enum AcquisitionWindowCalculationMode
%   {
%     AWCM_STANDARD            = 0x01,
%     AWCM_CONSIDER_LINES      = 0x02,
%     AWCM_CONSIDER_PARTITIONS = 0x04
%   };
% 
%   enum ExcitationPulse
%   {
%     EXCITATION_SLICE_SELECTIVE  = 0x01,
%     EXCITATION_VOLUME_SELECTIVE = 0x02
%   };
% 
%   enum SaturationRecovery
%   {
%     SATREC_NONE                  = 0x01,
%     SATREC_SLICE_SELECTIVE       = 0x02,
%     SATREC_VOLUME_SELECTIVE      = 0x04,
%     SATREC_VOLUME_SELECTIVE_PERF = 0x08,
%   };
% 
%   enum Device
%   {
%     DEVICE_DC  = 0x00,
%     DEVICE_RX  = 0x01,
%     DEVICE_TX  = 0x02,
%     DEVICE_ALL = 0xff
%   };
% 
%   enum ArrhythmiaDetection
%   {
%     AD_NONE         = 0x01,
%     AD_TIMEBASED    = 0x02,
%     AD_PATTERNBASED = 0x04
%   };
% 
%   enum FlowSensitivity
%   {
%     FLOW_SENSITIVITY_SLOW   = 0x01,
%     FLOW_SENSITIVITY_MEDIUM = 0x02,
%     FLOW_SENSITIVITY_FAST   = 0x04
%   };
% 
%   enum SharedDimension
%   {
%     SHARED_DIM_NONE         = 1,
%     SHARED_DIM_PHASES       = 2,
%     SHARED_DIM_SETS         = 3,
%     SHARED_DIM_REPETITIONS  = 4,
%     SHARED_DIM_ECHOES       = 5,
%     SHARED_DIM_FREE         = 6,
%     SHARED_DIM_ACQUISITIONS = 7,
%     SHARED_DIM_PARTITIONS   = 8
%   };
% 
%   enum ICEMode
%   {
%     ICE_MODE_MIP_SLAB_OVERLAP  = 0x01,
%     ICE_MODE_CISS              = 0x02,
%     ICE_MODE_DESS              = 0x04,
%     ICE_MODE_RE_DEPHASED       = 0x08
%   };
% 
%   enum FlowCompensation
%   {
%     FLOW_COMPENSATION_NO            = 0x01,
%     FLOW_COMPENSATION_YES           = 0x02,
%     FLOW_COMPENSATION_READOUT_ONLY  = 0x04,
%     FLOW_COMPENSATION_SLICESEL_ONLY = 0x08,
%     FLOW_COMPENSATION_SLICE_READ    = 0x10,
%     FLOW_COMPENSATION_SLICE_PHASE   = 0x20,
%     FLOW_COMPENSATION_READ_PHASE    = 0x40
%   };
% 
%   enum Tagging
%   {
%     TAGGING_NONE     = 0x01,
%     TAGGING_GRID_TAG = 0x02,
%     TAGGING_LINE_TAG = 0x04
%   };
% 
%   enum OnlineFFT
%   {
%     ONLINE_FFT_NONE      = 1,
%     ONLINE_FFT_PHASE     = 2,
%     ONLINE_FFT_PARTITION = 3
%   };
% 
%   enum FlowDir
%   {
%       FLOW_DIR_PHASE    = 0x01,
%       FLOW_DIR_READ     = 0x02,
%       FLOW_DIR_SLICESEL = 0x04,
%       FLOW_DIR_FREE     = 0x08,
%       FLOW_DIR_NONE     = 0x10        // not for protocol use!
%   };
% 
%   enum PCAngioFlowMode
%   {
%       PCANGIO_MODE_SINGLE_VELOCITY  = 0x01,
%       PCANGIO_MODE_SINGLE_DIRECTION = 0x02,
%       PCANGIO_MODE_FREE = 0x04
%   };
% 
%   enum PCAlgorithm
%   {
%     PC_ALGORITHM_NONE              = 1,
%     PC_ALGORITHM_SUBMATRIX         = 2,
%     PC_ALGORITHM_MARGOSIAN         = 3,
%     PC_ALGORITHM_3D_SUBMATRIX      = 4,
%     PC_ALGORITHM_3D_MARGOSIAN      = 5,
%     // the enum for POCS are defined in that way, that they also can be used as bits and masks
%     PC_ALGORITHM_POCS_RO           = 0x100,
%     PC_ALGORITHM_POCS_PE           = 0x200,
%     PC_ALGORITHM_POCS_3D           = 0x400,
%     PC_ALGORITHM_POCS_RO_PE        = 0x300,
%     PC_ALGORITHM_POCS_RO_3D        = 0x500,
%     PC_ALGORITHM_POCS_PE_3D        = 0x600,
%     PC_ALGORITHM_POCS              = 0x700
%   };
% 
%   enum NormalizeFilterAlgo
%   {
%     NORMALIZE_FILTER_ALGO_STANDARD       = 1,
%     NORMALIZE_FILTER_ALGO_HEAD           = 2,
%     NORMALIZE_FILTER_ALGO_ABDOMINAL      = 3
%   };
% 
%   enum ReorderingScheme
%   {
%     RS_LIN_IN_PAR          = 1,   // lines in partition
%     RS_PAR_IN_LIN          = 2,   // partitions in lines
%     RS_ARBITRARY           = 3    // arbitrary reordering schemes e.g. centric reordering
%   };
% 
%   enum PSatMode
%   {
%     PSAT_NONE         = 0x01,
%     PSAT_SINGLE_REG   = 0x02,
%     PSAT_DOUBLE_REG   = 0x04,
%     PSAT_SINGLE_QUICK = 0x08,
%     PSAT_DOUBLE_QUICK = 0x10
%   };
% 
%   enum RSatMode
%   {
%       RSAT_REG   = 0x01,
%       RSAT_QUICK = 0x02
%   };
% 
%   enum ShowOffline
%   {
%     SO_SHOW_YES          = 1,
%     SO_SHOW_NO           = 2
%   };
% 
%   enum FilterType
%   {
%     FILTER_NONE       = 0x01,
%     FILTER_RAW        = 0x02,        // RAW_FILTER clashes with prot.h
%     LARGE_FOV         = 0x04,
%     NORMALIZE         = 0x08,
%     ELLIPTICAL        = 0x10,
%     HAMMING           = 0x20,
%     FILTER_IMAGE      = 0x40,
%     PRESCAN_NORMALIZE = 0x80,
%     FILTER_BIFIC      = 0x100
%   };
% 
%   enum Reordering
%   {
%       REORDERING_LINEAR    = 0x01,
%       REORDERING_CENTRIC   = 0x02,
%       REORDERING_LINE_SEGM = 0x04,
%       REORDERING_PART_SEGM = 0x08,
%       REORDERING_FREE_0    = 0x10,
%       REORDERING_FREE_1    = 0x20,
%       REORDERING_FREE_2    = 0x40,
%       REORDERING_FREE_3    = 0x80
%   };
% 
%   enum PhaseStabScanPosition
%   {
%     AFTER        = 1,
%     BEFORE       = 2
%   };
% 
%   enum FlowDirDisplay
%   {
%     FLOW_DIR_R_TO_L     = R_TO_L,       // (R)ight (to) (L)eft
%     FLOW_DIR_L_TO_R     = L_TO_R,       // (L)eft (to) (R)ight
%     FLOW_DIR_A_TO_P     = A_TO_P,       // (A)nterior (to) (P)osterior
%     FLOW_DIR_P_TO_A     = P_TO_A,       // (P)osterior (to) (A)nterior
%     FLOW_DIR_F_TO_H     = F_TO_H,       // (F)eet (to) (H)ead
%     FLOW_DIR_H_TO_F     = H_TO_F,       // (H)ead (to) (F)eet
%     FLOW_DIR_TP_IN      = 0x020,        // (T)hrough (P)lane (In)flow
%     FLOW_DIR_TP_OUT     = 0x040,        // (T)hrough (P)lane (Out)flow
%     FLOW_DIR_INVALID    = 0x080         // No flow encoding, not for Protocol use
%   };
% 
% 
%   enum TriggerMode
%   {
%     TRIGGER_STANDARD     = 0x01,
%     TRIGGER_STEADY_STATE = 0x02
%   };
% 
%   enum ProtocolPackage
%   {
%     PROTPACKAGE_NONE      = 0x00,
%     PROTPACKAGE_MR_GUIDED = 0x01
%   };
% 
%   enum RfasSelMode
%   {
%     RFAS_SEL_NORMAL       = 0x00,
%     RFAS_SEL_TOGGLE_GAIN  = 0x01
%   };
% 
% 
%   //  Defines for dynamic image numbering
%   enum ImageNumbSag
%   {
%       IMAGE_NUMB_R_TO_L     = 0x0,
%       IMAGE_NUMB_L_TO_R     = 0x1,
%       IMAGE_NUMB_MED_TO_LAT = 0x2,
%       IMAGE_NUMB_LAT_TO_MED = 0x3
%   };
% 
%   enum ImageNumbCor
%   {
%       IMAGE_NUMB_A_TO_P = 0x0,
%       IMAGE_NUMB_P_TO_A = 0x1
%   };
% 
%   enum ImageNumbTra
%   {
%       IMAGE_NUMB_F_TO_H = 0x0,
%       IMAGE_NUMB_H_TO_F = 0x1
%   };
% 
%   enum ImageNumbMSMA
%   {
%       IMAGE_NUMB_SCT = 0x0,      // sagittal images precede coronal images,
%                                  // coronal images precede transverse images
%       IMAGE_NUMB_STC = 0x1,      // sagittal images precede transverse images,
%                                  // transverse images precede coronal images
%       IMAGE_NUMB_CTS = 0x2,      // coronal images precede transverse images,
%                                  // transverse images precede sagittal images
%       IMAGE_NUMB_CST = 0x3,      // coronal images precede sagittal images,
%                                  // sagittal images precede tranversal images
%       IMAGE_NUMB_TSC = 0x4,      // tranversal images precede sagittal images,
%                                  // sagittal images precede coronal images
%       IMAGE_NUMB_TCS = 0x5       // tranversal images precede coronal images,
%                                  // coronal images precede sagittal images
%   };
% 
%   enum DecouplingType
%   {
%     DECOUPLING_NONE    = 0x1,      // No decoupling
%     DECOUPLING_WALTZ4  = 0x2,      // WALTZ4 decoupling
%     DECOUPLING_MLEV    = 0x4,      // MLEV decoupling
%     DECOUPLING_CW      = 0x8,      // CW decoupling
%     DECOUPLING_WALTZ16 = 0x10,     // WALTZ16 decoupling
%   };
% 
%   enum NOEType
%   {
%     NOE_NONE         = 0x1,      // No NOE
%     NOE_RECTANGULAR  = 0x2       // Rectangular NOE shapes
%   };
% 
%   enum ExcitationType
%   {
%     EXCITATION_STANDARD  = 0x1,  // Standard excitation pulse type
%     EXCITATION_ADIABATIC = 0x2   // Adiabatic excitation pulse type
%   };
% 
%   enum PulseMode
%   {
%     EXCIT_MODE_2D_PENCIL            = 0x1,      // The timing used is a 2D-pencil-excitation
%     EXCIT_MODE_GRE                  = 0x2,      // The timing used is a gradient-echo sequence
%     EXCIT_MODE_EPI                  = 0x4,      // The timing uses a EPI sequence
%     EXCIT_MODE_CROSSED_PAIR         = 0x8,      // Two crossed slices with 90 and 180 degree excitations, forming a spin-echo where the slices intersect
%     EXCIT_MODE_2D_PENCIL_CARDIAC    = 0x10      // The timing used is a 2D-pencil-excitation for cardiac appls
% 
%   };
% 
%   enum ExamAnatomy
%   {
%     EXAM_ANATOMY_ABDOMEN = 0x01,    // selects an appropriate tracking factor for the abdomen
%     EXAM_ANATOMY_APEX    = 0x02,    // selects an appropriate tracking factor for the apex
%     EXAM_ANATOMY_LCA     = 0x03,    // selects an appropriate tracking factor for the LCA
%     EXAM_ANATOMY_OTHER   = 0x04     // the user can alter the tracking factor
%   };
% 
%   enum PATSelMode
%   {
%     PAT_MODE_NONE   = 0x01,
%     PAT_MODE_GRAPPA = 0x02,
%     PAT_MODE_SENSE  = 0x04,
%     PAT_MODE_2D     = 0x08
%   };
% 
%   enum PATRefScanMode
%   {
%     PAT_REF_SCAN_UNDEFINED      = 0x01, // e.g. if no PAT is selected
%     PAT_REF_SCAN_INPLACE        = 0x02, // sequence supplies inplace reference lines
%     PAT_REF_SCAN_EXTRA          = 0x04, // sequence supplies extra reference lines
%     PAT_REF_SCAN_PRESCAN        = 0x08, // sequence does not supply reference lines, the data must have been acquired with a previous measurement
%     PAT_REF_SCAN_INTRINSIC_AVE  = 0x10, // The sequence contains intrinsic ref.lines due to sharing e.g. in the averages dimension
%     PAT_REF_SCAN_INTRINSIC_REP  = 0x20, // The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
%     PAT_REF_SCAN_INTRINSIC_PHS  = 0x40, // The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
%     PAT_REF_SCAN_INPLACE_LET    = 0x80  // A single (L)ong (E)cho (T)rain acquires reference lines and imaging lines
%   };
% 
%   enum ChronPos
%   {
%     CHRON_POS_BEFORE_ECHO_TRAIN             = 0x01, // The navigator is actually executed before the image echo train
%     CHRON_POS_AFTER_ECHO_TRAIN              = 0x02, // The navigator is actually executed after the image echo train
%     CHRON_POS_BEFORE_AND_AFTER_ECHO_TRAIN   = 0x03  // The navigator is actually executed before and after the image echo train
%   };
% 
%   enum RspCompMode
%   {
%     RESP_COMP_BREATH_HOLD_AND_FOLLOW              =  0x01, // the navigator data is used to mitor the patient respiratory curve
%     RESP_COMP_GATE_AND_FOLLOW                     =  0x02, // navigator result is in a acceptance window, slices are shifted and measurement continues
%     RESP_COMP_OFF                                 =  0x04, // no respiratory compensation
%     RESP_COMP_GATE                                =  0x08, // data is aquired only if diaphragma position is in a acceptance window
%     RESP_COMP_BREATH_HOLD                         =  0x10, // slices are aquired after the operator pressed the scan button of the online display
%     RESP_COMP_BREATH_HOLD_AND_MONITOR             =  0x20, // similar to BREATH_HOLD_AND_FOLLOW, the differnce here is that the position from one breathhold to the next  is not adjusted
%     RESP_COMP_TRIGGER                             =  0x40, //
%     RESP_COMP_TRIGGER_AND_FOLLOW                  =  0x80, //
%     RESP_COMP_BREATH_HOLD_AND_TRIGGER             = 0x100, //
%     RESP_COMP_BREATH_HOLD_AND_TRIGGER_AND_FOLLOW  = 0x200, //
%     RESP_COMP_TRIGGER_AND_REACQ                   = 0x400, //
%     RESP_COMP_TRIGGER_AND_REACQ_AND_FOLLOW        = 0x800, //
%     RESP_COMP_MONITOR_ONLY                        = 0x1000 // The navigators are played out but they are not used for controlling the sequence
%   };
% 
%   enum BreathHoldNavigatorCapability
%   {
%     BREATHHOLD_NAVIGATOR_RERUN             = 0x01, // wether sequence allows user to rerun breathhold
%     BREATHHOLD_NAVIGATOR_SCAN              = 0x02, // wether sequence allows user to scan breathhold
%     BREATHHOLD_NAVIGATOR_DISPLAY_SELECTION = 0x04  // enables filter in Inline Dispaly, which either blocks secondary captured images or non-secondary captured images
%   };
% 
%   enum DiffusionDirectionality
%   {
%     DIFFDIR_NONE         = 0x01,    // none specifies diffusion conditions
%     DIFFDIR_DIRECTIONAL  = 0x02,    // specifies whether diffusion conditions for the frame are directional with respect to direction
%     DIFFDIR_ISOTROPIC    = 0x04     // specifies whether diffusion conditions for the frame are isotropic with respect to direction
%   };
% 
%   enum TRFillPosition
%   {
%     BEFORE_ACQUISITION = 1,     // the delay time for each acquisition (TRFill) is before the acquisition
%     AFTER_ACQUISITION  = 2      // the delay time for each acquisition (TRFill) is after the acquisition
%   };
% 
%   enum PALIMode
%   {
%     PALI_MODE_NORMAL     = 1,       // normal PALI mode for patient examinations
%     PALI_MODE_SERVICE    = 2        // PALI disabled (for some service sequences only, not for patient examinations!!!!)
%   };
% 
%   //---------------------------------------------------------------------------
%   // Allows the sequence to control details of the application functionality - see
%   // parameter ulApplicationDetails in prot.h
%   //
%   // - the concrete functionality usually depends on ulApplicationDetails plus the
%   //   selected postprocessing protocol (parameter tDefaultEVAProt in prot.h)
%   //---------------------------------------------------------------------------
%   enum ApplicationDetails
%   {
%     APPLICATION_NONE          = 0x00, //deactivates all special application parameter settings
%     APPLICATION_INLINE_BREAST = 0x01, //activates the breast perfusion card in case the inline+breast postprocessing is active
%     APPLICATION_IRT           = 0x02, //controls I(interactive R(ealT(ime dependencies to XProtocol based parameters
%     APPLICATION_TIMCT         = 0x04  //controls TIM C(ontinue T(able move dependencies to XProtocol based parameters
%   };
% 
%   enum ApplSpec
%   {
%       APPL_BODY      = 0x01,//Application for Spectroscopy
%       APPL_PROSTATE  = 0x02  
%   };
% 
%   enum SpectralSuppression
%   {
%     SPEC_SUPPR_NONE         = 0x01,  // no spectral suppression mode
%     SPEC_SUPPR_LIPID        = 0x02,  // lipid spectral suppression mode
%     SPEC_SUPPR_WATER        = 0x04,  // water spectral suppression mode
%     SPEC_SUPPR_LIPID_WATER  = 0x06   // lipid AND water spectral suppression modes
%   };
% 
%   enum SegmentationMode
%   {
%     SEGM_MODE_DEFINE_SEGM_FACTOR        = 0x01,  //  
%     SEGM_MODE_DEFINE_SHOTS              = 0x02   //  
%   };
% 
%   enum Trajectory
%   {
%     TRAJECTORY_CARTESIAN  = 0x01,
%     TRAJECTORY_RADIAL     = 0x02,
%     TRAJECTORY_SPIRAL     = 0x04,
%     TRAJECTORY_BLADE      = 0x08,
%   };
% 
%   enum ViewSharing
%   {
%     VIEW_SHARING_OFF     = 0x01,
%     VIEW_SHARING_SHPHS   = 0x02,
%     VIEW_SHARING_KTBLAST = 0x04,
%     VIEW_SHARING_TWIST   = 0x08
%   };
% 
%   enum SequenceType
%   {
%     SEQUENCE_TYPE_UNDEF  = 0x00,
%     SEQUENCE_TYPE_GRE    = 0x01,
%     SEQUENCE_TYPE_TRUFI  = 0x02,
%     SEQUENCE_TYPE_EPI    = 0x04,
%     SEQUENCE_TYPE_TSE    = 0x08,
%     SEQUENCE_TYPE_CSI    = 0x10,
%     SEQUENCE_TYPE_FID    = 0x20,
%   };
% 
%   enum CoilCombineMode
%   {
%     COILCOMBINE_SUM_OF_SQUARES      = 0x01,
%     COILCOMBINE_ADAPTIVE_COMBINE    = 0x02,
%   };
% 
%   enum FlipAngleMode
%   {
%     FAM_CONSTANT           = 0x01,
%     FAM_HYPERECHO          = 0x02,
%     FAM_SAR_OPTIMIZED_T1   = 0x04,
%     FAM_SAR_OPTIMIZED_PD   = 0x08,
%     FAM_SAR_OPTIMIZED_T2   = 0x10,
%     FAM_VARIABLE           = 0x20,
%   };
%       
%   enum TOM
%   {
%     TOM_OFF                  = 0x01,
%     TOM_MINIMIZE_TE          = 0x02,
%     TOM_MINIMIZE_TR          = 0x04,
%     TOM_MINIMIZE_TE_TR       = 0x08,
%     TOM_MINIMIZE_ECHOSPACING = 0x10,
%     TOM_IN_PHASE             = 0x20,
%     TOM_OPPOSED_PHASE        = 0x40,
%   };
% 
%   enum MultipleSeries
%   {
%     MULTIPLE_SERIES_OFF                        = 0x01,
%     MULTIPLE_SERIES_EACH_SLICE                 = 0x02,  // images of every slice are put into a new series
%     MULTIPLE_SERIES_EACH_MEASUREMENT           = 0x04,  // images of every measurement/repetition are put into a different series
%     MULTIPLE_SERIES_EACH_SLICE_AND_MEASUREMENT = 0x08,  // create a new series for each slice and repetition.  
%   };
% 
% 
%   enum DistortionCorrMode
%   {
%     DISTCORR_NDIS           = 0x01,
%     DISTCORR_DIS2D          = 0x02,
%     DISTCORR_DIS3D          = 0x04,
%   };
%   
%   enum TablePositioningMode
%   {
%     TP_POS_MODE_FIX         = 0x01,
%     TP_POS_MODE_GSP         = 0x02,
%     TP_POS_MODE_ISO         = 0x04,
%   };
%   
%   enum ReadOutMode
%   {
%       READOUT_MONOPOLAR     = 0x01,
%       READOUT_BIPOLAR       = 0x02,
%   };
% 
%   enum AsymmetricEchoMode
%   {
%       ASYMM_ECHO_WEAK       = 0x1,
%       ASYMM_ECHO_STRONG     = 0x2,
%       ASYMM_ECHO_HALF       = 0x4
%   };
% 
%   enum BladeMotionCorr
%   {
%       BLADE_MOTION_CORR_ON          = 0x01, // creates a motion-corrected series
%       BLADE_MOTION_CORR_OFF         = 0x02, // creates a not-motion-corrected series
%       BLADE_MOTION_CORR_ON_OFF      = 0x04, // creates two series, one motion corrected an one not corrected
%   };
% 
%   enum MDSMode
%   {
%        MDS_OFF   = 0x01,
%        MDS_ON    = 0x02,
%        MDS_ONCO  = 0x04
%   };
% 
%   enum MDSVariableResolution
%   {
%        MDS_VARIABLE_RESOLUTION_OFF          = 0x01,
%        MDS_VARIABLE_RESOLUTION_END_OF_RANGE = 0x02
%   };
% 
%   enum MDSReconMode
%   {
%        MDS_RECON_TOTAL              = 0x01,
%        MDS_RECON_CHOPPED            = 0x02,
%        MDS_RECON_TOTAL_AND_CHOPPED  = 0x04
%   };
% 
%   enum AdjShim
%   {
%     ADJSHIM_TUNEUP      = 0x01,    // use tuneup shim setting
%     ADJSHIM_STANDARD    = 0x02,    // perform standard shim procedure
%     ADJSHIM_ADVANCED    = 0x04,    // perform advanced shim procedure
%     ADJSHIM_INTERACTIVE = 0x20,    // perform interactive shim procedure (internal use)
%     ADJSHIM_MEAS_ONLY   = 0x80,    // measurement only, no calculation (internal use)
%     ADJSHIM_CALC_ONLY   = 0x40     // calculation only, no measurement (internal use)
%   };
% 
%   enum AdjFre
%   {
%     ADJFRE_FINE           = 0x01,  // use narrow band adj/fre for exact results
%     ADJFRE_COARSE         = 0x02,  // use wide band adj/fre for coarse results
%     ADJFRE_FINE_NOTVOLSEL = 0x04   // use narrow band adj/fre for exact results (not volume selective)
%   };
% 
%   enum AdjMDS
%   {
%     ADJMDS_NONE   = 0x01,  // no MDS activities
%     ADJMDS_ADJUST = 0x02,  // measure MDS adjustments
%     ADJMDS_SCOUT  = 0x04,  // measure MDS scout
%     ADJMDS_BOTH   = 0x06,  // measure MDS adjustments and scout
%   };
% 
%   enum AdjustWithTolerance
%   {
%     ADJTOL_NONE             = 0x01,
%     ADJTOL_MAXIMUM          = 0x02,
%   };
%   
%   enum ParametricMapMode
%   {
%        PMAP_NONE        = 0x01,
%        PMAP_T1_MAP      = 0x02,
%        PMAP_T2_MAP      = 0x04,
%        PMAP_T2STAR_MAP  = 0x08
%   };
% 
%   enum AngioDynamicReconMode
%   {
%       ANGIO_DYNAMIC_RECON_FORWARD_SHARE                 = 0x01,
%       ANGIO_DYNAMIC_RECON_BACKWARD_SHARE                = 0x02,
%       ANGIO_DYNAMIC_RECON_LINEAR_INTERPOLATION          = 0x04,
%       ANGIO_DYNAMIC_RECON_CUBIC_SPLINE_INTERPOLATION    = 0x08
%   };
%  
%   ///Selection of ASL mode to measure perfusion.
%   enum AslMode
%   {
%     ASL_NONE            = 0x01,
%     ASL_PICOREQ2TIPS    = 0x02,
%     ASL_FAIRQUIPSSII    = 0x04,
%     ASL_PICOREQUIPSSII  = 0x08,
%     ASL_CASL            = 0x10,
%     ASL_PSEUDOCASL      = 0x20,
%     ASL_CUSTOM1         = 0x40,
%     ASL_CUSTOM2         = 0x80
% 
%   };
%   
%   /// tx excitation mode
%   enum BCExcitationMode
%   {
%         /// invalid/empty value
%         TX_BC_INVALID         = 0x00,
%         /// always use CP mode, regardless of coils
%         TX_BC_CP              = 0x01,
%         /// always use elliptic mode, regardless of coils
%         TX_BC_ELLIPTIC        = 0x02,
%         /// use auto mode (determined by coil select)
%         TX_BC_AUTO            = 0x04
%   };
% 
%   enum TmapB0Correction
%   {
%       TMAP_B0_CORRECTION_ON  = 0x01,
%       TMAP_B0_CORRECTION_OFF = 0x02
%   };
% 
%   enum TmapEval
%   {
%       TMAP_EVAL_ON  = 0x01,
%       TMAP_EVAL_OFF = 0x02  
%   };
% 
%   enum TmapImageType
%   {
%       TMAP_REF_IMAGE        = 0x01,
%       TMAP_TEMPERTURE_IMAGE = 0x02
%   };
% 
%   //---------------------------------------------------------------------------
%   // Defines can be used in two ways:
%   //
%   // 1. SEQ::NUCLEUS_1H  -> SEQ::MeasNucleus("1H")
%   //                     -> SEQ_MEASNUCLEUS("1H")
%   //                     -> MeasNucleus("1H")
%   // 2. NUCLEUS_1H       -> MeasNucleus("1H")
%   //---------------------------------------------------------------------------
%   #ifdef __cplusplus
%     typedef SEQ_MEASNUCLEUS     MeasNucleus;
%     typedef SEQ_MEASNUCLEUSSET  MeasNucleusSet;
%     typedef SEQ_MEASKNOWNNUCLEI MeasKnownNuclei;
%   #endif
% 
%   #define NUCLEUS_1H          MeasNucleus("1H")
%   #define NUCLEUS_3HE         MeasNucleus("3He")
%   #define NUCLEUS_7LI         MeasNucleus("7Li")
%   #define NUCLEUS_13C         MeasNucleus("13C")
%   #define NUCLEUS_17O         MeasNucleus("17O")
%   #define NUCLEUS_19F         MeasNucleus("19F")
%   #define NUCLEUS_23NA        MeasNucleus("23Na")
%   #define NUCLEUS_31P         MeasNucleus("31P")
%   #define NUCLEUS_129XE       MeasNucleus("129Xe")
%   #define NUCLEUS_NONE        MeasNucleus("")        // Nucleus not defined
% 
%   #define NUCLEI_NONE         MeasNucleusSet("")     // Empty set
%   #define NUCLEI_ALL          MeasKnownNuclei()      // All known nuclei
% };
% 
% 
% //-----------------------------------------------------------------------------
% // Definitions for array lengths (used in SeqLim.h)
% //-----------------------------------------------------------------------------
% #define MAX_SEQ_LIM_FILENAMELENGTH    128
% #define MAX_SEQ_LIM_INT_PAR           16
% #define MAX_SEQ_LIM_FLOAT_PAR         16
% 
% //-----------------------------------------------------------------------------
% // Definitions for array lengths (used in SeqExpo.h, YAPS)
% //-----------------------------------------------------------------------------
% #define MAX_SUPPORTED_MAXWELL_COEFFICIENTS  64
% 
% //------------------------------------------------------------------------------
% // Nametags for Breathhold Navigator used by Sequence and Online Display
% //------------------------------------------------------------------------------
% const char BREATHOLD_NAMETAG_SCAN[]  =   "BH_SCAN";
% const char BREATHOLD_NAMETAG_RERUN[] =   "BH_RRUN";
% #endif
% 
% /*---------------------------------------------------------------------------*/
% /*  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential        */
% /*---------------------------------------------------------------------------*/

% /*---------------------------------------------------------------------------*/
% /*  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential        */
% /*---------------------------------------------------------------------------*/
% /*
%  * Project: NUMARIS/4
%  *    File: \n4\pkg\MrServers\MrMeasSrv\SeqIF\MDH\mdh.h
%  * Version:
%  *  Author: CC_MEAS SCHOSTZF
%  *    Date: n.a.
%  *
%  *    Lang: C
%  *
%  * Descrip: measurement data header
%  *
%  *                      ATTENTION
%  *                      ---------
%  *
%  *  If you change the measurement data header, you have to take care that
%  *  long variables start at an address which is aligned for longs. If you add
%  *  a short variable, then add two shorts from the short array or use the
%  *  second one from the previous change (called "dummy", if only one was added and
%  *  no dummy exists).
%  *  Then you have to extend the swap-method from MdhProxy.
%  *  This is necessary, because a 32 bit swaped is executed from MPCU to image
%  *  calculator.
%  *  Additional, you have to change the dump method from libIDUMP/IDUMPRXUInstr.cpp.
%  *
%  * Functns: n.a.
%  *
%  *---------------------------------------------------------------------------*/
% 
% /*--------------------------------------------------------------------------*/
% /* Include control                                                          */
% /*--------------------------------------------------------------------------*/
% #ifndef MDH_H
% #define MDH_H
% 
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
% #define MDH_NUMBEROFICEPROGRAMPARA 4
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of free header parameters (short)                            */
% /*--------------------------------------------------------------------------*/
% #define MDH_FREEHDRPARA  (4)
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
% } sLoopCounter;                                /* sizeof : 28 byte             */
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of slice vectors                                             */
% /*--------------------------------------------------------------------------*/
% 
% typedef struct
% {
%   PACKED_MEMBER( float,  flSag          );
%   PACKED_MEMBER( float,  flCor          );
%   PACKED_MEMBER( float,  flTra          );
% } sVector;
% 
% typedef struct
% {
%   sVector sSlicePosVec; /* slice position vector        */
%   PACKED_MEMBER( float,           aflQuaternion[4] ); /* rotation matrix as quaternion*/
% } sSliceData;                                         /* sizeof : 28 byte             */
% 
% /*--------------------------------------------------------------------------*/
% /*  Definition of cut-off data                                              */
% /*--------------------------------------------------------------------------*/
% typedef struct
% {
%   PACKED_MEMBER( uint16_t,  ushPre          );    /* write ushPre zeros at line start */
%   PACKED_MEMBER( uint16_t,  ushPost         );    /* write ushPost zeros at line end  */
% } sCutOffData;
% 
% 
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
% #include "MrServers/MrVista/include/unpack.h"
% 
% #endif   /* MDH_H */
% 
% /*---------------------------------------------------------------------------*/
% /*  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential        */
% /*---------------------------------------------------------------------------*/



%% ---------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------

%% VD13



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

tic
S0 = size(Data);
maxCOL = S0(1);
numOfLines = S0(2);

maxCOL = 0;
startInd = double(asc(1).ushUsedChannels+1);
for ii=startInd:numOfLines    
    if ( asc(ii).ushSamplesInScan > maxCOL )
        maxCOL = asc(ii).ushSamplesInScan;
    end
end

% find the noise lines and seperate ref lines if any
numOfNoiseLines = 0; 
NoiseLines = zeros(numOfLines,1, 'int32');
sLCNoise = zeros(numOfLines, 11, 'int32');

numOfSeperateRefLines = 0; 
RefLines = zeros(numOfLines,1, 'int32');
numOfRefAndImgLines = 0; 
RefAndImgLines = zeros(numOfLines,1, 'int32');
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

numOfOtherLines = 0; 
OtherLines = zeros(numOfLines,1, 'int32');
sLCOther = zeros(numOfLines, 11, 'int32');

kspaceCentreLineNo = zeros(numOfLines,1, 'int32');

currREPForOther = 1;
for i=1:numOfLines
      
    if ( IsAcqEnd(asc(i).aulEvalInfoMask(1)) )
        continue;
    end
    
    kspaceCentreLineNo(i) = asc(i).ushKSpaceCentreLineNo;
    
    [bNoiseLine, bSeperateRef, bRefandImg, bPhaseCorr, bReflect, bLastScanInSlice] = parseEvalInfoMask(asc(i).aulEvalInfoMask(1));
    
    if ( bNoiseLine  )
        numOfNoiseLines = numOfNoiseLines + 1;
        NoiseLines(numOfNoiseLines) = i;
        sLCNoise(numOfNoiseLines,:) = fillSLC(asc(i));
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
    
    if ( bSeperateRef || bRefandImg )
        numOfSeperateRefLines = numOfSeperateRefLines + 1;
        RefLines(numOfSeperateRefLines) = i;
        sLCRef(numOfSeperateRefLines,:) = fillSLC(asc(i));
        
        if ( bSeperateRef )
            continue;
        end
        
        if ( bRefandImg )
            numOfRefAndImgLines = numOfRefAndImgLines + 1;
            RefAndImgLines(numOfRefAndImgLines) = i;
        end
    end
        
    if ( asc(i).ushSamplesInScan < maxCOL & ~bRefandImg  )
        numOfOtherLines = numOfOtherLines + 1;
        OtherLines(numOfOtherLines) = i;
        sLCOther(numOfOtherLines,:) = fillSLC(asc(i));
        sLCOther(numOfOtherLines,9) = currREPForOther;
        continue;
    end
        
    if ( bReflect ) % if it is a reflected data lines
        numOfReflectLines = numOfReflectLines + 1;
        ReflectLines(numOfReflectLines) = i;
        sLCReflect(numOfReflectLines,:) = fillSLC(asc(i));
    end 

    numOfDataLines = numOfDataLines + 1; % include normal data line and refAndImg lines
    DataLines(numOfDataLines) = i;
    sLC(numOfDataLines,:) = fillSLC(asc(i));

    if ( bLastScanInSlice & sLC(numOfDataLines,5)==1 )
        currREPForOther = sLC(numOfDataLines,9) + 1;
        % disp(['current REP : ' num2str(currREPForOther)]);
    end
end

sLCNoise = sLCNoise(1:numOfNoiseLines, :); NoiseLines = NoiseLines(1:numOfNoiseLines);

sLCRef = sLCRef(1:numOfSeperateRefLines, :); RefLines = RefLines(1:numOfSeperateRefLines); RefAndImgLines = RefAndImgLines(1:numOfRefAndImgLines);

sLCPhaseCorr = sLCPhaseCorr(1:numOfPhaseCorrLines, :); PhaseCorrLines = PhaseCorrLines(1:numOfPhaseCorrLines);

sLCReflect = sLCReflect(1:numOfReflectLines, :); ReflectLines = ReflectLines(1:numOfReflectLines);

sLC = sLC(1:numOfDataLines, :); DataLines = DataLines(1:numOfDataLines);

sLCOther = sLCOther(1:numOfOtherLines, :); OtherLines = OtherLines(1:numOfOtherLines);

maxKSpaceLineNo = 2*max(kspaceCentreLineNo);
[kspace, kspaceASCIndex] = allocateKSpace(sLC, maxKSpaceLineNo);

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

% other lines
if ( numOfOtherLines > 0 )
    [other, otherASCIndex] = allocateKSpace(sLCOther, maxKSpaceLineNo);
else
    other = [];
    otherASCIndex = [];
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
    
    % fill other lines
    r = find(i==OtherLines);
    if ( ~isempty(r)  )
        N = asc(i).ushSamplesInScan;
        aZ = addPrePostZeros(asc(i));
        r = r(1);
        if ( sLCOther(r,2) <= maxKSpaceLineNo )
                       
            if ( aZ == 0 )
                other(1:N, sLCOther(r,2), sLCOther(r,3), sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i);
            elseif ( aZ == 1 ) % pre zeros
                other(end-N+1:end, sLCOther(r,2), sLCOther(r,3), sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i);
            elseif ( aZ == 2 ) % post zeros
                other(1:N, sLCOther(r,2), sLCOther(r,3), sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = Data(1:N, i);
            end
            
            otherASCIndex(sLCOther(r,2), sLCOther(r,3), sLCOther(r,4), sLCOther(r,5), ...
                    sLCOther(r,6), sLCOther(r,7), sLCOther(r,8), sLCOther(r,9), sLCOther(r,10), sLCOther(r,11)) = i;
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
                        
            if ( aZ == 0 )
                ref(1:N, sLCRef(r,2), sLCRef(r,3), sLCRef(r,4), sLCRef(r,5), ...
                    sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = Data(1:N, i);
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
        
        pr = find(i==RefAndImgLines);
        if ( isempty(pr)  )
            continue;
        end
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
    function [bNoiseLine, bSeperateRef, bRefAndImg, bPhaseCorr, bReflect, bLastScanInSlice] = parseEvalInfoMask(evalInfoMask)
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

