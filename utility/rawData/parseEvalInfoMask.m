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
    