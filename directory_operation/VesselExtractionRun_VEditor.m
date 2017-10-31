
function [vasculature, handles] = VesselExtractionRun_VEditor(data, header, presetscales, handles, seeds, scales)
% perform the vessel extraction

% extracting vessels

% global definition and declaration
global vasculature
global Aylward_coefficient
% vasculature
vasculature = struct('count',0, 'vessels', [], 'totalnumofpoints', 0);
onevessel = struct('ridgepoints',[], 'radius',[], 'samplepositions', [], 'len', 0, 'dist', [],...
    'points_eigenvalue', [], 'points_eigenvector1', [], 'points_eigenvector2', [], ...
    'points_eigenvector3', [], 'points_sigma', []);
vasculatureNodes = struct('count',0, 'nodes', [], 'totalnumofpoints', 0);

samplestep = 0.2;
numofkernels = 5;
%-------------------------------------------------------------------------%
doradius = 0;

samplestepnew = 0.05;
presetScales = [0.2 : 0.2 : 4];


% compute the coefficient
if ( isempty(handles.coefficient) == 1 )
    handles.coefficient = col2row( ComputeCoefficients(data) );
end

if ( handles.otherparameters.DynamicScale ~= 1.0 ) % use default scale
%     vasculature = RidgeDetection5_SeedsScales(data, header, samplestepnew, handles.parameters, ...
%         handles.radiusparameters, doradius, handles.coefficient, seeds, scales, handles.coefficientBlurred, handles.otherparameters.DefaultScale);
    
    vasculature = RidgeDetection5_SeedsScales_noPerturbation(data, header, samplestepnew, handles.parameters, ...
        handles.radiusparameters, doradius, handles.coefficient, seeds, scales, handles.coefficientBlurred, handles.otherparameters.DefaultScale);

else
%     vasculature = RidgeDetection5_SeedsScales(data, header, samplestepnew, handles.parameters, ...
%         handles.radiusparameters, doradius, handles.coefficient, seeds, scales, handles.coefficientBlurred, -1);
    
    vasculature = RidgeDetection5_SeedsScales_noPerturbation(data, header, samplestepnew, handles.parameters, ...
        handles.radiusparameters, doradius, handles.coefficient, seeds, scales, handles.coefficientBlurred, -1);

end

% save vasculature into the extractionResults directory...
vesselfile = fullfile('.', 'vasculature_VEditor_seeds.mat');
save(vesselfile, 'vasculature');
return;