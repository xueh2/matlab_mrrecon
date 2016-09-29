function kspace = PerformRampSamplingOffCenterCorrection(kspace, asc, kspaceASCIndex, reflect, headers, protocol)
% for the ramp sampling, perform the off center correction

if ( isempty(kspace) | isempty(reflect) )
    return;
end

% find the rampup, rampdown, flattop, delaytime, readoutFOV
[feFOV, peFOV, sliceThickness] = findFOVFromConfig(headers.Config);
[RampupTime, RampdownTime, FlattopTime, DelaySamplesTime, ADCDuration, DestSamples] = findRampSamplingTimes(headers.MrcTime);

if ( RampupTime == 0 )
    disp(['No ramp sampling ... ']);
    return;
end

RampupTime_in_s = RampupTime*1e-6; % in second
RampdownTime_in_s = RampdownTime*1e-6;
FlattopTime_in_s = FlattopTime*1e-6;
DelaySamplesTime_in_s = DelaySamplesTime*1e-6;
ADCDuration_in_s = ADCDuration*1e-6;

COL = size(kspace, 1);
% rx_dwelltime = protocol.sRXSPEC.alDwellTime_{1};

% deltaT_in_s = 1/(10^9/COL/rx_dwelltime);
deltaT_in_s = ADCDuration_in_s/COL;

% compute the idea kspace sampling location
t = 0:COL;
t = t * deltaT_in_s + DelaySamplesTime_in_s;

k_sampled = zeros(COL+1, 1);
for i=1:COL+1
    if ( t(i) < RampupTime_in_s )
        k_sampled(i) = 0.5*t(i)*t(i)/RampupTime_in_s;
    elseif ( (t(i)>=RampupTime_in_s) & (t(i)<RampupTime_in_s+FlattopTime_in_s) )
        k_sampled(i) = 0.5*RampupTime_in_s + t(i)-RampupTime_in_s;
    else
        dt = t(i) - RampupTime_in_s - FlattopTime_in_s;
        gt = (RampdownTime_in_s - dt)/RampdownTime_in_s;
        k_sampled(i) = 0.5*RampupTime_in_s + FlattopTime_in_s + 0.5*dt*(gt+1);
    end
end    

% scaleFac = COL/(k_sampled(end)-k_sampled(1));
% k_sampled = k_sampled * scaleFac;

% compute the acquired kspace sampling location, to make sure the kspace center matches
k_linear = k_sampled;
for i=COL/2:-1:1
    k_linear(i) = k_linear(i+1) - deltaT_in_s; 
end

for i=COL/2+2:COL+1
    k_linear(i) = k_linear(i-1) + deltaT_in_s; 
end

% k_linear = 0.5*DelaySamplesTime_in_s*DelaySamplesTime_in_s/RampupTime_in_s + (t-DelaySamplesTime_in_s);
% k_linear = reshape(k_linear, [COL 1]);

% k_linear = 1:COL;

% scale the sampling location, so the sinc interpolation makes sense
scalefac = COL / (k_sampled(end)-k_sampled(1));%( k_linear(end) - k_linear(1) ) / ( k_sampled(end) - k_sampled(1) );
k_sampled = scalefac * k_sampled;
k_linear = COL / (k_linear(end)-k_linear(1)) * k_linear;

% hold on; plot(k_sampled, 'r.'); plot(k_linear, 'b.'); hold off
 
clear i

% compute the phase offset
[COL, LIN, CHA, ACQ, SLC, PAR, ECO, PHS, REP, SET, SEG] = findKSpaceDimension(kspace);

fReadOutOffcentrePre = -1;
phaseOffset = zeros(COL, 1);
for seg=1:SEG
    for set=1:SET
        for rep=1:REP
            for phs=1:PHS
                for eco=1:ECO
                    for par=1:PAR
                        for slc=1:SLC                        
                            for acq=1:ACQ
                                for lin=1:LIN

                                    ascInd = kspaceASCIndex(lin, 1, acq, slc, par, eco, phs, rep, set, seg);
                                    if ( ascInd > 0 )
                                        fReadOutOffcentre = asc(ascInd).fReadOutOffcentre;

                                        if ( fReadOutOffcentre ~= fReadOutOffcentrePre )
                                            % need to compute phase offset
                                            rFOV = fReadOutOffcentre/(2*feFOV); % 2x readout oversampling
                                            phaseOffset = 2*pi*rFOV.*(k_linear(1:COL) - k_sampled(1:COL));
                                            phaseOffset = exp(-2i*pi.*phaseOffset);
                                            phaseOffset = reshape(phaseOffset, [COL 1]);
                                            phaseOffset = repmat(phaseOffset, [1 CHA]);
                                            phaseOffset = reshape(phaseOffset, [COL 1 CHA]);
                                        end

                                        reflectFlag = reflect(1, lin, 1, acq, slc, par, eco, phs, rep, set, seg);
                                        if ( reflectFlag )
                                            readOut = kspace(:,lin,:, acq, slc, par, eco, phs, rep, set, seg);
                                            readOut = flipdim(readOut, 1);
                                            readOut = readOut .* phaseOffset;
                                            kspace(:,lin,:, acq, slc, par, eco, phs, rep, set, seg) = flipdim(readOut, 1);
                                        else
                                            kspace(:,lin,:, acq, slc, par, eco, phs, rep, set, seg) = kspace(:,lin,:, acq, slc, par, eco, phs, rep, set, seg) .* phaseOffset;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
