function out = excitation_profile(bt, plotflag);
% function excitation_profile(bt);
%
% function to plot the slice/slab excitation profile for a Hanning
% weighted sinc as a function of timebandwidth product
%
% Usage:
%     excitation_profile([3.5]); % single bt value
%     excitation_profile([1.6 2 4 6 8]); % multiple bt values

if nargin < 2; plotflag = 1; end

if plotflag; hfig = figure; hold on; end

for i=1:length(bt)
    N=1000; % number of samples per RF pulse
    timebandwidthproduct = bt(i);
    x = linspace(-timebandwidthproduct/2, timebandwidthproduct/2, N*timebandwidthproduct);
    rf = sinc(x).*hanning(length(x))'/N;
    M=200; % zerofilled interpolation factor
    profile(:,i) = fftshift(fft(rf, M*length(rf)/timebandwidthproduct));
end
f = linspace(-N/2, N/2, size(profile,1));
if plotflag==1;
    plot(f, abs(profile));
    % draw prescribed slice (vertical dotted lines)
    plot([-.5 -.5],[ 0 1.2],'k:')
    plot([.5 .5],[ 0 1.2],'k:')
    % set axes
    axis([-4 4 0 1.2])
    box on
    legend (num2str([bt(:)]));
    shg
elseif plotflag==2
    plot(f, abs(profile),'w');
    % draw prescribed slice (vertical dotted lines)
    plot([-.5 -.5],[ 0 1.2],'w:')
    plot([.5 .5],[ 0 1.2],'w:')
    % set axes
    axis([-1 1 0 1.2])
    box on
    bwaxis
    h_legend=legend (num2str([bt(:)]));
    set(h_legend,'color','k','edgecolor','w','textcolor','w')
    shg
elseif plotflag == 3
    plot(f, abs(profile));
    hold on;
    factor = 1/1.398; % empirical factor
    plot(factor * f, abs(profile), 'r');
    % draw prescribed slice (vertical dotted lines)
    plot([-.5 -.5],[ 0 1.2],'k:')
    plot([.5 .5],[ 0 1.2],'k:')
    % set axes
    axis([-1 1 0 1.2])
    box on
    legend (num2str([bt(:)]));
    shg    
end

if nargout > 0
    index_min = max(find(f <= -2.5));
    index_max = min(find(f >=  2.5));
    out.profile = abs(profile(index_min:index_max, :));
    out.z = f(index_min:index_max);
end
    



    
    
% % default 3d
% m_RFProp.getDuration (ePulseType)  800
% m_RFProp.getBWT      (ePulseType)  3.5
% m_RFProp.getEmpiricalFactor()      1
% 
% % default 2d
% m_RFProp.getDuration (ePulseType)  600
% m_RFProp.getBWT      (ePulseType)  1.6
% m_RFProp.getEmpiricalFactor()      1


% excerpt from a_trufi_cv.cpp
% // * -------------------------------------------------------------------------- *
% // *                                                                            *
% // * Name        :  Trufi_cv::p_SelectRFPulse                                   *
% // *                                                                            *
% // * Description :  Parametrization of the RF pulses                            *
% // *                This function is new in a_CV                                *
% // *                                                                            *                                                                            *
% // * Return      :  NLS status                                                  *
% // *                                                                            *
% // * -------------------------------------------------------------------------- *
% NLS_STATUS Trufi_cv::p_SelectRFPulse ( 
%                                      MrProt           *pMrProt,       // * IMP: Measurement protocol *
%                                      SeqLim           *pSeqLim,       // * IMP: Sequence limits      *
%                                      SeqExpo          *               // * IMP: Returned values      *
%                                      )
% {
%    
%     static const char *ptModule = { "Trufi_cv::p_SelectRFPulse" };
%     mTRrun;
% 
%                                  // (Duration FNL   ,     BWT-Prod. FNL,     empFactor)             
%     // cine trueFISP/turboFlash
%     RFProperties RFCine2DTFL        (1000, 2000, 4000,   2.0 , 2.7,  2.0,    1.0      ); // based on a_gre_retro
%     RFProperties RFCine2DTFLTrio    (1000, 2000, 4000,   2.0 , 2.7,  2.0,    1.0      ); // based on a_gre_retro
% 
%     RFProperties RFCine2DTFI        ( 600, 1000, 1600,   1.6 , 1.6,  1.6,    1.0      ); // based on cine_trufi_retro
%     RFProperties RFCine2DTFIPre     ( 600, 1000, 1600,   1.6 , 1.6,  1.6,    1.0      ); // based on cine_trufi_retro
% 
%     RFProperties RFCine2DTFITrio    ( 800, 1200, 1600,   1.6 , 1.6,  1.6,    0.8      ); // based on cine_trufi_retro
%     RFProperties RFCine2DTFITrioPre ( 800,  800,  800,   1.6 , 1.6,  1.6,    0.8      ); // based on cine_trufi_retro
% 
%     // dynamic imaging only available for 2D sequences, no prep pulse for tfl mode
%     RFProperties RFTFLPerf2D        ( 400,  800, 1200,   1.6 , 2.0,  2.6,    1.0      ); // based on a_tfl, low SAR not available in a_tfl 
% 
%     RFProperties RFTFIPerf2D        ( 600, 1000, 1400,   1.6 , 2.0,  2.6,    1.0      ); // based on a_tfiperf, low SAR not available in a_tfiperf
%     RFProperties RFTFIPerf2DPre     ( 500,  800, 1200,   1.6 , 2.0,  2.6,    1.0      ); // based on a_tfiperf, different length due to PNS
%     
%     // magprep 
%     RFProperties RFTFLSeg2D         ( 600, 1000, 1400,   1.6 , 2.0,  2.6,    1.0      ); // based on     a_irtfl2D 
% 
%     //RFProperties RFTFISeg2D         ( 600, 1000, 1400,   1.6 , 2.0,  2.6,    1.0      ); // based on  a_tfiseg, low SAR not available in a_tfiseg;
%     RFProperties RFTFISeg2D         ( 600, 1200, 1400,   1.6 , 4.0,  2.6,    1.0      ); // PK*
%     RFProperties RFTFISeg2DPre      ( 500,  800, 1200,   1.6 , 2.0,  2.6,    1.0      ); // based on  a_tfiseg, different length due to PNS
% 
% 
%     // 3D pulses vvvvvv
%     // cine trueFISP/turboFlash 
%     RFProperties RFStd3D            (1000, 1500, 2000,   6.0 , 6.4,  6.4,    1.0      ); // UNUSED! just for reference
% 
%     RFProperties RFCine3DTFL        ( 600,  800, 1400,   6.0,  9.0, 12.0,    1.0      ); 
%     RFProperties RFCine3DTFLTrio    ( 600,  800, 1400,   6.0,  9.0, 12.0,    1.0      ); 
% 
%     RFProperties RFCine3DTFI        ( 600, 1000, 1600,   6.0,  9.0, 12.0,    1.0      ); 
%     RFProperties RFCine3DTFIPre     ( 600, 1000, 1600,   6.0,  9.0, 12.0,    1.0      ); 
% 
%     RFProperties RFCine3DTFITrio    ( 600, 1000, 1600,   6.0 , 9.0, 12.0,    1.0      ); 
%     RFProperties RFCine3DTFITrioPre ( 600, 1000, 1600,   6.0 , 9.0, 12.0,    1.0      ); 
%  
%     // dynamic imaging
%     RFProperties RFTFLPerf3D        ( 600, 1000, 1600,   3.5 , 6.4, 10.2,    1.0      ); // based on a_tfiseg
% 
%     RFProperties RFTFIPerf3D        ( 600, 1000, 1600,   3.5 , 6.4, 10.2,    1.0      ); // based on a_tfiperf, low SAR not available in a_tfiperf
%     RFProperties RFTFIPerf3DPre     ( 500,  800, 1200,   3.5 , 6.4,  9.6,    1.0      ); // based on a_tfiperf, different length due to PNS
% 
%     // magprep
%     RFProperties RFTFLSeg3D         (1000, 2000, 4000,   6.4 ,12.7, 25.4,    1.0      ); // based on a_tfl,     low SAR not available in a_tfl 
% 
%     RFProperties RFTFISeg3D         ( 800, 1000, 1400,   3.5 , 4.5,  6.4,    1.0      ); // based on a_tfiseg, Charm 346452     
%     RFProperties RFTFISeg3DPre      ( 800, 1000, 1400,   3.5 , 4.5,  6.4,    1.0      ); // based on a_tfiseg, length due to PNS
% 
%     RFProperties RFTFISeg3DTrio     (1000, 1200, 1600,   3.5 , 4.5,  6.4,    1.0      ); // based on a_tfiseg  Charm 346452  
%     RFProperties RFTFISeg3DTrioPre  (1000, 1200, 1600,   3.5 , 4.5,  6.4,    1.0      ); // based on a_tfiseg, length due to PNS
% 
%     RFProperties RFnon_sel          ( 250,  300,  400,   1.0 , 1.0,  1.0,    1.0      ); // non-selective excitation  
%     RFProperties RFnon_selTrio      ( 400,  500,  600,   1.0 , 1.0,  1.0,    1.0      ); // non-selective excitation
% 
%     RFProperties RFnon_selUTE       (  60,  100,  400,   1.0 , 1.0,  1.0,    1.0      ); // non-selective excitation
% 

