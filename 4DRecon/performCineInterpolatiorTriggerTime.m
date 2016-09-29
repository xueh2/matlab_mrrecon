function phaseIms = performCineInterpolatiorTriggerTime(data, triggerTime, cardiacPhases, method)
% function phaseIms = performCineInterpolatiorTriggerTime(data, triggerTime, cardiacPhases)
%
% this function performs the retrospective gating by interpolating the cardiac phase images at specific trigger time
%
% Inputs:
%    data : [Nfe Npe numel(triggerTime)], the input cardiac phase images
%    triggerTime : trigger times of input cardiac phase images
%    cardiacPhases : requried cardiac phases
%
% Output:
%    phaseIms : interpolated cardiac phase images
%
%     ***************************************
%     *  Hui Xue (hui-xue@siemens.com       *
%     *  2012-03                            *
%     ***************************************
if ( nargin<4 )
    method = 'linear';
end

if ( strcmp(method, 'linear') )
    Nfe = size(data, 1);
    Npe = size(data, 2);
    NPhs = numel(cardiacPhases);
    inputNPhs = numel(triggerTime);

    phaseIms = zeros(Nfe, Npe, NPhs);

    for k=1:NPhs

        phs = cardiacPhases(k);
        if ( phs <= triggerTime(1) )
            phaseIms(:,:,k) = data(:,:,1);
            continue;
        end

        if ( phs >= triggerTime(end) )
            phaseIms(:,:,k) = data(:,:,end);
            continue;
        end

        for t=1:inputNPhs-1
            if ( phs>=triggerTime(t) & phs<triggerTime(t+1)   )
                break;
            end
        end

        alpha = (phs-triggerTime(t))/(triggerTime(t+1)-triggerTime(t));
        phaseIms(:,:,k) = (1-alpha)*data(:,:,t) + alpha*data(:,:,t+1);
    end
end

if ( strcmp(method, 'pchip') )
    Nfe = size(data, 1);
    Npe = size(data, 2);
    NPhs = numel(cardiacPhases);
    inputNPhs = numel(triggerTime);

    phaseIms = zeros(Nfe, Npe, NPhs);

    for pe=1:Npe
        for fe=1:Nfe
            phaseIms(fe,pe,:) = interp1(triggerTime', squeeze(data(fe,pe,:)), cardiacPhases', 'pchip');
        end
    end

    for k=1:NPhs    
        phs = cardiacPhases(k);
        if ( phs <= triggerTime(1) )
            phaseIms(:,:,k) = data(:,:,1);
            continue;
        end

        if ( phs >= triggerTime(end) )
            phaseIms(:,:,k) = data(:,:,end);
            continue;
        end   
    end
end

% %% apply the FFD based key frame estimation
% 
% pixelSpacing = [1 1 1];
% header = CreateFtkHeaderInfo(data, [1 1 1]);
% N = size(data, 3)
% realVolume = cell(N, 1);
% imagVolume = cell(N, 1);
% headers = cell(N, 1);
% for i=1:N    
%     realVolume{i} = double(real(data(:,:,i)));
%     imagVolume{i} = double(imag(data(:,:,i)));
%     header2D = Dicom2HeaderMrFtk(realVolume{i}, pixelSpacing, [0 0 triggerTime(i)], [1 0 0 0 1 0]);
%     if ( i ~= N ) 
%         header2D.spacingZ = triggerTime(i+1)-triggerTime(i);
%     else
%         header2D.spacingZ = 1;
%     end
%     headers{i} = header2D;
% end
% 
% % ensure the maximal level of approximation
% numOfSubdivision = 9;
% controlPointSpacing = 0.5*[header.spacingX*header.sizeX  header.spacingY*header.sizeY  triggerTime(end)-triggerTime(1)];
% 
% % recon every phase
% NPhs = numel(cardiacPhases);
% 
% headerDst = header;
% headerDst.sizeZ = NPhs+4;
% headerDst.spacingZ = (cardiacPhases(end)-cardiacPhases(1))/(NPhs-1);
% headerDst.positionPatient = [0 0 cardiacPhases(1)-2*headerDst.spacingZ];
% headerDst.orientationPatient = eye(3,3);
% 
% for i=1:headerDst.sizeZ
%     [wx, wy, wz ] = Image2WorldMrFtk(headerDst, 0, 0, i-1);
%     wz
% end
% 
% realPhaseIms = Matlab_ScatterInterpolationFFD(realVolume, headers, headerDst, controlPointSpacing, numOfSubdivision, []);
% imagPhaseIms = Matlab_ScatterInterpolationFFD(imagVolume, headers, headerDst, controlPointSpacing, numOfSubdivision, []);
% 
% phaseIms = complex(realPhaseIms, imagPhaseIms);
