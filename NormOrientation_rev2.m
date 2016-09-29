function  [rotate90, flip] = NormOrientation_rev2 (rowvec, colvec, normvec);
% function  [rotate90, flip] = NormOrientation_rev2 (quaternion);
% 
% main orientation is SAG, COR, or TRA
%
% operations to perform: flip then rotate90
% rotate90 = 0, +/-1 (for use in matlab rot90.m function)
%            1 (CCW 90 == CW 270) , -1 (CW 90)
% flip = 0  (no flip)
%      = 1  (flip rows only)
%      = 2  (flip cols only)

%                                 standard            	   swapped
% sagittal                        rot 90 (CW)                 -
% coronal                         rot 270, mirr COL       mirr LIN
% transversal, cranial (default)  mirr COL                rot 270, mirr COL
% transversal, caudal             rot 270

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************


verbose=0;% for debugging

% R = quat2mat(quaternion); % convert quaternion to rotation matrix
% 
% colvec  = R(1,:); % 1st row of rotation matrix
% rowvec  = R(2,:); % 2nd row of rotation matrix
% normvec = R(3,:); % 3rd row of rotation matrix

% determine image main orientations of normal vector & row vector
[orientation, value] = MainOrientation (rowvec, colvec, normvec);

mainElementRow = value.row;
mainElementColumn = value.col;
mainOrientation=orientation.main;
rowOrientation= orientation.row;

%   // Basis of row and column vector is the patient coordinate system (SAG, COR, TRA).
%   // Intput is an image coming with any orthogonal row and column vector.
%   // This image shall be oriented in such way, that
%   // * sagittal   images are shown with row vector being positive and column vector being negative (+-),
%   // * coronal    images are shown with row vector being positive and column vector being negative (+-),
%   // * transverse images are shown with row vector being positive and column vector being positive (++).
%   // A vector being 'positive' or 'negative' means that it's biggest element is positive or negative.
%   // The main orientation is determined by the biggest element of the slice normal vector.
%   // Depending on the orientation (++ or +-) different operations have to be performed for desired result.

true=1; false=0;
isOrientedRow = false;
isOrientedColumn = false;

switch(mainOrientation)
    case 1 % SAG
        desiredOrientationRow = 2;
    case 2 % COR
        desiredOrientationRow = 1;
    case 3 % TRA
        desiredOrientationRow = 1;
end

switch(mainOrientation)
    case {1,2} %SAG, COR
    %     // -----------------------------
    %     // - case for orientation (+-) -
    %     // -----------------------------
        if(rowOrientation == desiredOrientationRow)
    %         // -------------------
    %         // - 'parallel' case -
    %         // -------------------
            isOrientedRow = mainElementRow>0;%       // row compared to row direction (+)
            isOrientedColumn = mainElementColumn<0;% // column compared to columnn direction (-)
        else
    %         // ---------------------
    %         // - 'orthogonal' case -
    %         // ---------------------
            isOrientedRow = mainElementRow<0; %      // row compared to column direction (-)
            isOrientedColumn = mainElementColumn>0;% // column compared to row direction (+)
        end
    case 3 % TRA
%               // -----------------------------
%               // - case for orientation (++) -
%               // -----------------------------
            isOrientedRow = mainElementRow>0;
            isOrientedColumn = mainElementColumn>0;
end
              
% // control rotate/mirror of images
rotationAngle = 0;
mirrorDimension = -1;
      
if(rowOrientation == desiredOrientationRow)
% // -------------------
% // - 'parallel' case -
% // -------------------
    if verbose==1
        disp(['parallel case:', num2str(rowOrientation)])
    end
    if(isOrientedRow && isOrientedColumn)
    end
    if(isOrientedRow && ~isOrientedColumn)
        mirrorDimension = 1; % // mirror LIN
    end
    if(~isOrientedRow && isOrientedColumn)
        mirrorDimension = 0; %// mirror COL
    end
    if(~isOrientedRow && ~isOrientedColumn)
        rotationAngle = 180;
%         rotationAngle = 0; % hack
    end
else
%   // ---------------------
%   // - 'orthogonal' case -
%   // ---------------------
    if verbose==1
        disp(['orthogonal case:', num2str(rowOrientation)])
    end
    if(mainOrientation==1) % SAG
    % // In this case we have to invert our flags to have identical cases in if clauses.
    % // sag: (cor x -tra is antiparallel to sag) <- !!
    % // cor: (sag x -tra is     parallel to cor)
    % // tra: (sag x  cor is     parallel to tra)
        isOrientedRow = ~isOrientedRow;
        isOrientedColumn = ~isOrientedColumn;
    end
    % // invert flags again if imageNormal points to positive main direction
    % // (original images from measurement always have negative main direction for normal)
    % // used to avoid reverse rotation
    normvec(1) = rowvec(2)*colvec(3)-rowvec(3)*colvec(2);
    normvec(2) = rowvec(3)*colvec(1)-rowvec(1)*colvec(3);
    normvec(3) = rowvec(1)*colvec(2)-rowvec(2)*colvec(1);
    index=find(abs(normvec)==max(abs(normvec)));
    MainElementNormVec = normvec(index(1));

    if(MainElementNormVec>0)
        isOrientedRow = ~isOrientedRow;
        isOrientedColumn = ~isOrientedColumn;
    end
    if(isOrientedRow && isOrientedColumn)
        rotationAngle = 90;
        mirrorDimension = 1; %// mirror LIN
    end
    if(~isOrientedRow && isOrientedColumn)
        rotationAngle = 270;
    end
    if(isOrientedRow && ~isOrientedColumn)
        rotationAngle = 90;
    end
    if(~isOrientedRow && ~isOrientedColumn)
        rotationAngle = 270;
        mirrorDimension = 1; % // mirror LIN
    end
end

switch rotationAngle
    case 0
        rotate90 = 0; % do nothing
    case 90
        rotate90 = 1;% 90 deg CCW
    case 180
        rotate90 = 2; % 180
    case 270
        rotate90 = -1; % 90 deg CW or 270 CCW
end
switch mirrorDimension
    case -1
        flip = 0; % do nothing
    case 0
        flip = 2; % mirror COL   -> flip left/right
    case 1
        flip = 1; % mirror LIN   -> flip up/down    
end

return



% from IceSodaFactory.h
%         double rotMatrix[3][3];
%         double* Gp = rotMatrix[0]; // The GP vector
%         double* Gr = rotMatrix[1]; // The GR vector
%         double* Gs = rotMatrix[2]; // The GS vector (= slice normal vector)
% 
%         dstSoda._rowVec.dSag = Gr[0];
%         dstSoda._rowVec.dCor = Gr[1];
%         dstSoda._rowVec.dTra = Gr[2];
% 
%         dstSoda._columnVec.dSag = Gp[0];
%         dstSoda._columnVec.dCor = Gp[1];
%         dstSoda._columnVec.dTra = Gp[2];
        
% IceWrapper2.cpp
% bool IceWrapper::NormOrientation

%   // Basis of row and column vector is the patient coordinate system (SAG, COR, TRA).
%   // Intput is an image coming with any orthogonal row and column vector.
%   // This image shall be oriented in such way, that
%   // * sagittal   images are shown with row vector being positive and column vector being negative (+-),
%   // * coronal    images are shown with row vector being positive and column vector being negative (+-),
%   // * transverse images are shown with row vector being positive and column vector being positive (++).
%   // A vector being 'positive' or 'negative' means that it's biggest element is positive or negative.
%   // The main orientation is determined by the biggest element of the slice normal vector.
%   // Depending on the orientation (++ or +-) different operations have to be performed for desired result.

%           switch(mainOrientation)
%           {
%           case IceSoda::SAG:
%               desiredOrientationRow = IceSoda::COR;
%               break;
%           case IceSoda::COR:
%               desiredOrientationRow = IceSoda::SAG;
%               break;
%           case IceSoda::TRA:
%               desiredOrientationRow = IceSoda::SAG;
%               break;
%           case IceSoda::NONE:
%               ICE_ERR("NormOrientation was called with invalid SliceOrientationData."
%                       "Could not determine the mainOrientation.");
%               ICE_ERR("Normalize Orientation failed because rotate/mirror reported an error.");
%               IceFramework::SetLastError(I_FAIL);
%               return (false);
%           }
%           // -------------------------------
%           // - switch over mainOrientation -
%           // -------------------------------
% 
%           bool isOrientedRow = false;
%           bool isOrientedColumn = false;
% 
%           switch(mainOrientation)
%           {
%           case IceSoda::SAG:
%           case IceSoda::COR:
% 
%               // -----------------------------
%               // - case for orientation (+-) -
%               // -----------------------------
% 
%               if(rowOrientation == desiredOrientationRow)
%               {
%                   // -------------------
%                   // - 'parallel' case -
%                   // -------------------
%                   isOrientedRow = mainElementRow>0;       // row compared to row direction (+)
%                   isOrientedColumn = mainElementColumn<0; // column compared to columnn direction (-)
%               }
%               else
%               {
%                   // ---------------------
%                   // - 'orthogonal' case -
%                   // ---------------------
%                   isOrientedRow = mainElementRow<0;       // row compared to column direction (-)
%                   isOrientedColumn = mainElementColumn>0; // column compared to row direction (+)
%               }
% 
%               break;
%           case IceSoda::TRA:
% 
%               // -----------------------------
%               // - case for orientation (++) -
%               // -----------------------------
% 
%               isOrientedRow = mainElementRow>0;
%               isOrientedColumn = mainElementColumn>0;
% 
%               break;
%           }
% 
%           if(rowOrientation == desiredOrientationRow)
%           {
%               // -------------------
%               // - 'parallel' case -
%               // -------------------
% 
%               if(isOrientedRow && isOrientedColumn)
%               {
%                   // do nothing
%               }
% 
%               if(isOrientedRow && !isOrientedColumn)
%               {
%                   mirrorDimension = 1; // mirror LIN
%               }
% 
%               if(!isOrientedRow && isOrientedColumn)
%               {
%                   mirrorDimension = 0; // mirror COL
%               }
% 
%               if(!isOrientedRow && !isOrientedColumn)
%               {
%                   rotataionAngle = SDCRotation180;
%               }
%           }
%           else
%           {
%               // ---------------------
%               // - 'orthogonal' case -
%               // ---------------------
%               if(mainOrientation==IceSoda::SAG)
%               {
%                   // In this case we have to invert our flags to have identical cases in if clauses.
%                   // sag: (cor x -tra is antiparallel to sag) <- !!
%                   // cor: (sag x -tra is     parallel to cor)
%                   // tra: (sag x  cor is     parallel to tra)
%                   isOrientedRow = !isOrientedRow;
%                   isOrientedColumn = !isOrientedColumn;
%               }
% 
%               // invert flags again if imageNormal points to positive main direction
%               // (original images from measurement always have negative main direction for normal)
%               // used to avoid reverse rotation
%               sVec imageNormal;
%               imageNormal.dSag = rowVec.dCor*columnVec.dTra-rowVec.dTra*columnVec.dCor;
%               imageNormal.dCor = rowVec.dTra*columnVec.dSag-rowVec.dSag*columnVec.dTra;
%               imageNormal.dTra = rowVec.dSag*columnVec.dCor-rowVec.dCor*columnVec.dSag;
%               if(imageNormal.getMainElement()>0)
%               {
%                   isOrientedRow = !isOrientedRow;
%                   isOrientedColumn = !isOrientedColumn;
%               }
% 
%               if(isOrientedRow && isOrientedColumn)
%               {
%                   rotataionAngle = SDCRotation90;
%                   mirrorDimension = 1; // mirror LIN
%               }
% 
%               if(!isOrientedRow && isOrientedColumn)
%               {
%                   rotataionAngle = SDCRotation270;
%               }
% 
%               if(isOrientedRow && !isOrientedColumn)
%               {
%                   rotataionAngle = SDCRotation90;
%               }
% 
%               if(!isOrientedRow && !isOrientedColumn)
%               {
%                   rotataionAngle = SDCRotation270;
%                   mirrorDimension = 1; // mirror LIN
%               }
%           }
%       }



