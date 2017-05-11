function  [orientation, value] = MainOrientation_rotation_matrix(ImageOrientationPatient);
% function  [orientation, value] = MainOrientation_rotation_matrix (ImageOrientationPatient);
%
% quaternion
%
% main orientation is SAG, COR, or TRA

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

% R = quat2mat(quaternion); % convert quaternion to rotation matrix

% colvec  = R(1,:); % 1st row of rotation matrix
% rowvec  = R(2,:); % 2nd row of rotation matrix
% normvec = R(3,:); % 3rd row of rotation matrix

colvec  = ImageOrientationPatient(1:3); % 1st row of rotation matrix
rowvec  = ImageOrientationPatient(4:6); % 2nd row of rotation matrix
normvec = cross(rowvec,colvec); % 3rd row of rotation matrixs


% determine main image orientation from normal vector
% find direction with max direction
[dMaxNorm, maxOriNorm] = getMaxElement(normvec);
orientation.main = maxOriNorm;
value.main = dMaxNorm;

% index=find(abs(normvec)==max(abs(normvec)));
% switch index
%     case 1
%         orientation.main = 'SAG';
%     case 2
%         orientation.main = 'COR';
%     case 3
%         orientation.main = 'TRA';
% end
% // fallback, should be overwritten in every case!
dMainOriRow = rowvec(1);
dMainOriCol = colvec(1);

[dMinRow, minOriRow] = getMinElement(rowvec);
[dMaxRow, maxOriRow] = getMaxElement(rowvec);
[dMinCol, minOriCol] = getMinElement(colvec);
[dMaxCol, maxOriCol] = getMaxElement(colvec);

dDeltaRow = 0; dDeltaCol = 0; % // delta of 2 *largest* elements
% compute delta of the 'big ones'
switch(minOriRow)
    case 3 % TRA
        dDeltaRow = abs(abs(rowvec(1))-abs(rowvec(2)));
    case 2 % COR
        dDeltaRow = abs(abs(rowvec(1))-abs(rowvec(3)));
    case 1 % SAG
        dDeltaRow = abs(abs(rowvec(3))-abs(rowvec(2)));
end
switch(minOriCol)
    case 3 % TRA
        dDeltaCol = abs(abs(colvec(1))-abs(colvec(2)));
    case 2 % COR
        dDeltaCol = abs(abs(colvec(1))-abs(colvec(3)));
    case 1 % SAG
        dDeltaCol = abs(abs(colvec(3))-abs(colvec(2)));
end

% define true & false
true=1;
false=0;

% used to determine main orientation for small deltas
bRowSagIsMain = true;
bRowCorIsMain = true;
bRowTraIsMain = true;
bColSagIsMain = true;
bColCorIsMain = true;
bColTraIsMain = true;

% row and columns' orientation cannot be overall main orientation
switch(orientation.main)
    case 3 % TRA
        bRowTraIsMain = false;
        bColTraIsMain = false;
    case 2 % COR
        bRowCorIsMain = false;
        bColCorIsMain = false;
    case 1 % SAG
        bRowSagIsMain = false;
        bColSagIsMain = false;
end

% determine main orientations
% trust vector with larger delta
if(dDeltaRow>dDeltaCol) %main orientation of column vec cannot be main orientation of row
    switch(maxOriRow)
        case 2 % COR
            bColCorIsMain = false;
            dMainOriRow = rowvec(2);
            mainOriRow = 2;
        case 1 % SAG
            bColSagIsMain = false;
            dMainOriRow = rowvec(1);
            mainOriRow = 1;
        case 3 % TRA
            bColTraIsMain = false;
            dMainOriRow = rowvec(3);
            mainOriRow = 3;
    end
else %main orientation of row vec cannot be main orientation of column
    switch(maxOriCol)
        case 2 % COR
            bRowCorIsMain = false;
            dMainOriCol = colvec(2);
        case 1 % SAG
            bRowSagIsMain = false;
            dMainOriCol = colvec(1);
        case 3 % TRA
            bRowTraIsMain = false;
            dMainOriCol = colvec(3);
    end
end

if(dDeltaRow>dDeltaCol) % determine main element of column
    if(bColSagIsMain)
        dMainOriCol = colvec(1);
    end
    if(bColCorIsMain)
        dMainOriCol = colvec(2);
    end
    if(bColTraIsMain)
        dMainOriCol = colvec(3);
    end
else % determine main element of row
    if(bRowSagIsMain)
        dMainOriRow = rowvec(1);
        mainOriRow = 1;
    end
    if(bRowCorIsMain)
        dMainOriRow = rowvec(2);
        mainOriRow = 2;
    end
    if(bRowTraIsMain)
        dMainOriRow = rowvec(3);
        mainOriRow = 3;
    end
end

% output row and col orientations
orientation.row = mainOriRow;
value.row = dMainOriRow;
value.col = dMainOriCol;



function [dMinElement, minOrientation] = getMinElement(vector);
% function that outputs value of element in vector which has min magnitude (abs)
% minOrientation = 1,2,3 (Sag, Cor, Tra)
% index=find(abs(vector)==min(abs(vector)));
% value=vector(index(1));
% orientation=index;

if(abs(vector(3))<abs(vector(2)))
    dMinElement = vector(3);
	minOrientation = 3;
	if(abs(vector(1))<abs(vector(3)))
        dMinElement = vector(1);
        minOrientation = 1;
    end
else
    dMinElement = vector(2);
    minOrientation = 2;
    if(abs(vector(1))<abs(vector(2)))
        dMinElement = vector(1);
        minOrientation= 1;
    end
end


function [dMaxElement, maxOrientation] = getMaxElement(vector);
% function that outputs value of element in vector which has max magnitude (abs)
% maxOrientation = 1,2,3 (Sag, Cor, Tra)

% index=find(abs(vector)==max(abs(vector)));
% value=vector(index(1));
% orientation=index;

if(abs(vector(3))>abs(vector(2)))
    dMaxElement = vector(3);
	maxOrientation = 3;
	if(abs(vector(1))>abs(vector(3)))
        dMaxElement = vector(1);
        maxOrientation = 1;
    end
else
    dMaxElement = vector(2);
    maxOrientation = 2;
    if(abs(vector(1))>abs(vector(2)))
        dMaxElement = vector(1);
        maxOrientation= 1;
    end
end



return


  
% void IceSDC::getMainElements(const sVec& rowVec, const sVec& colVec, double& dMainOriRow, IceSoda::PatientKosy& mainOriRow, double& dMainOriCol, IceSoda::PatientKosy mainOrientation)
% {
%     ICE_SET_FN("IceSDC::getMainElements(const sVec& rowVec, const sVec& colVec, double& dMainOriRow, IceSoda::PatientKosy& mainOriRow, double& dMainOriCol, IceSoda::PatientKosy mainOrientation)")
% 
%     // We know the mainOrientation of the image (determined by classOri using the normal vector).
%     // We have to determine the main element of row and column vector.
%     // Note that getMainOrientation (i.e. classOri) must not be called using row/column vectors, it only works fine for normal vectors!
%     // We know that the main orientation of the image is neither main element of row nor column vector!
%     // In case of identity of 2 elements of a vector we try to use information retrieved from the other vector.
% 
%     // retrieve delta of largest elements of a vector and orientation of minimum
%     double dMinRow, dMaxRow;
%     IceSoda::PatientKosy minOriRow, maxOriRow;
%     double dMinCol, dMaxCol;
%     IceSoda::PatientKosy minOriCol, maxOriCol;
% 
%     // fallback, should be overwritten in every case!
%     dMainOriRow = rowVec.dSag;
%     dMainOriCol = colVec.dSag;
% 
%     // get min/max values and orientations
%     getMinElement(rowVec, dMinRow, minOriRow);
%     getMinElement(colVec, dMinCol, minOriCol);
%     getMaxElement(rowVec, dMaxRow, maxOriRow);
%     getMaxElement(colVec, dMaxCol, maxOriCol);
% 
%     double dDeltaRow = 0, dDeltaCol = 0; // delta of 2 *largest* elements
% 
%     // compute delta of the 'big ones'
%     switch(minOriRow)
%     {
%     case IceSoda::TRA:
%         dDeltaRow = fabs(fabs(rowVec.dSag)-fabs(rowVec.dCor));
%         break;
%     case IceSoda::COR:
%         dDeltaRow = fabs(fabs(rowVec.dSag)-fabs(rowVec.dTra));
%         break;
%     case IceSoda::SAG:
%         dDeltaRow = fabs(fabs(rowVec.dTra)-fabs(rowVec.dCor));
%         break;
%     default:
%         ICE_ERR("IceSDC::getMainElements: invalid minOriRow given (" << minOriRow << "), proceeding...");
%     }
% 
%     switch(minOriCol)
%     {
%     case IceSoda::TRA:
%         dDeltaCol = fabs(fabs(colVec.dSag)-fabs(colVec.dCor));
%         break;
%     case IceSoda::COR:
%         dDeltaCol = fabs(fabs(colVec.dSag)-fabs(colVec.dTra));
%         break;
%     case IceSoda::SAG:
%         dDeltaCol = fabs(fabs(colVec.dTra)-fabs(colVec.dCor));
%         break;
%     default:
%         ICE_ERR("IceSDC::getMainElements: invalid minOriCol given (" << minOriCol << "), proceeding...");
%     }
% 
%      // used to determine main orientation for small deltas
%     bool bRowSagIsMain = true;
%     bool bRowCorIsMain = true;
%     bool bRowTraIsMain = true;
%     bool bColSagIsMain = true;
%     bool bColCorIsMain = true;
%     bool bColTraIsMain = true;
% 
%     // row and columns' orientation cannot be overall main orientation
%     switch(mainOrientation)
%     {
%     case IceSoda::TRA:
%         bRowTraIsMain = false;
%         bColTraIsMain = false;
%         break;
%     case IceSoda::COR:
%         bRowCorIsMain = false;
%         bColCorIsMain = false;
%         break;
%     case IceSoda::SAG:
%         bRowSagIsMain = false;
%         bColSagIsMain = false;
%         break;
%     default:
%         ICE_ERR("IceSDC::getMainElements: invalid mainOrientation given (" << mainOrientation << "), proceeding...");
%     }
% 
%     // determine main orientations
%     // trust vector with larger delta
%     if(dDeltaRow>dDeltaCol)
%     {
%         // main orientation of column vec cannot be main orientation of row
%         switch(maxOriRow)
%         {
%         case IceSoda::COR:
%             bColCorIsMain = false;
%             dMainOriRow = rowVec.dCor;
%             mainOriRow = IceSoda::COR;
%             break;
%         case IceSoda::SAG:
%             bColSagIsMain = false;
%             dMainOriRow = rowVec.dSag;
%             mainOriRow = IceSoda::SAG;
%             break;
%         case IceSoda::TRA:
%             bColTraIsMain = false;
%             dMainOriRow = rowVec.dTra;
%             mainOriRow = IceSoda::TRA;
%             break;
%         default:
%             ICE_ERR("IceSDC::getMainElements: invalid maxOriRow (" << maxOriRow << "), proceeding...");
%         }
%     }
%     else
%     {
%         // main orientation of row vec cannot be main orientation of column
%         switch(maxOriCol)
%         {
%         case IceSoda::COR:
%             bRowCorIsMain = false;
%             dMainOriCol = colVec.dCor;
%             break;
%         case IceSoda::SAG:
%             bRowSagIsMain = false;
%             dMainOriCol = colVec.dSag;
%             break;
%         case IceSoda::TRA:
%             bRowTraIsMain = false;
%             dMainOriCol = colVec.dTra;
%             break;
%         default:
%             ICE_ERR("IceSDC::getMainElements: invalid maxOriCol (" << maxOriCol << "), proceeding...");
%         }
%     }
% 
%     if(dDeltaRow>dDeltaCol)
%     {
%         // determine main element of column
%         if(bColSagIsMain)
%         {
%             dMainOriCol = colVec.dSag;
%         }
%         if(bColCorIsMain)
%         {
%             dMainOriCol = colVec.dCor;
%         }
%         if(bColTraIsMain)
%         {
%             dMainOriCol = colVec.dTra;
%         }
%     }
%     else
%     {
%         // determine main element of row
%         if(bRowSagIsMain)
%         {
%             dMainOriRow = rowVec.dSag;
%             mainOriRow = IceSoda::SAG;
%         }
%         if(bRowCorIsMain)
%         {
%             dMainOriRow = rowVec.dCor;
%             mainOriRow = IceSoda::COR;
%         }
%         if(bRowTraIsMain)
%         {
%             dMainOriRow = rowVec.dTra;
%             mainOriRow = IceSoda::TRA;
%         }
%     }
% }
% 
% 
% 
% 
% 
% 
% void IceSDC::getMinElement(const sVec& vec, double& dMinElement, IceSoda::PatientKosy& minOrientation)
% {
%     ICE_SET_FN("IceSDC::getMinElement(const sVec& vec, double& dMinElement, IceSoda::PatientKosy& minOrientation)")
% 
%     if(fabs(vec.dTra)<fabs(vec.dCor))
%     {
%         dMinElement = vec.dTra;
%         minOrientation = IceSoda::TRA;
% 
%         if(fabs(vec.dSag)<fabs(vec.dTra))
%         {
%             dMinElement = vec.dSag;
%             minOrientation = IceSoda::SAG;
%         }
%     }
%     else
%     {
%         dMinElement = vec.dCor;
%         minOrientation = IceSoda::COR;
% 
%         if(fabs(vec.dSag)<fabs(vec.dCor))
%         {
%             dMinElement = vec.dSag;
%             minOrientation= IceSoda::SAG;
%         }
%     }
% }
% 
% 
% 
% void IceSDC::getMaxElement(const sVec& vec, double& dMaxElement, IceSoda::PatientKosy& maxOrientation)
% {
%     ICE_SET_FN("IceSDC::getMaxElement(const sVec& vec, double& dMaxElement, IceSoda::PatientKosy& maxOrientation)")
% 
%     if(fabs(vec.dTra)>fabs(vec.dCor))
%     {
%         dMaxElement = vec.dTra;
%         maxOrientation = IceSoda::TRA;
% 
%         if(fabs(vec.dSag)>fabs(vec.dTra))
%         {
%             dMaxElement = vec.dSag;
%             maxOrientation = IceSoda::SAG;
%         }
%     }
%     else
%     {
%         dMaxElement = vec.dCor;
%         maxOrientation = IceSoda::COR;
% 
%         if(fabs(vec.dSag)>fabs(vec.dCor))
%         {
%             dMaxElement = vec.dSag;
%             maxOrientation= IceSoda::SAG;
%         }
%     }
% }