function [dataCorr, posCorr, rowCorr, colCorr] = performFrameOrientationCorrection(data, pixelSpacing, pos, rowVector,colVector, posPre, rowVectorPre,colVectorPre)
% correct the frame orientation and roate/flip the data if needed

% create header
header = CreateFtkHeaderInfo(data, pixelSpacing);
header.positionPatient = pos;
header.orientationPatient(1,:) = rowVector;
header.orientationPatient(2,:) = colVector;
header.orientationPatient(3,:) = cross(header.orientationPatient(1,:), header.orientationPatient(2,:));

headerPre = CreateFtkHeaderInfo(data, pixelSpacing);
headerPre.positionPatient = posPre;
headerPre.orientationPatient(1,:) = rowVectorPre;
headerPre.orientationPatient(2,:) = colVectorPre;
headerPre.orientationPatient(3,:) = cross(headerPre.orientationPatient(1,:), headerPre.orientationPatient(2,:));

row = size(data, 1);
col = size(data, 2);

% compute the four corner of input images
[wc1, wc2, wc3, wc4] = computeFourCorner(data, header);

% compute posCorr
dist = [ norm(wc1-posPre); norm(wc2-posPre); norm(wc3-posPre); norm(wc4-posPre)];

[minDist, ind] = min(dist);
if ind==1
    posCorr = wc1;
    
    rowCorr = wc2-wc1;
    colCorr = wc4-wc1;
    
    cost = dot(rowVectorPre, rowCorr) + dot(colVectorPre, colCorr);
    if ( cost < 0.1 )
        temp = rowCorr;
        rowCorr = colCorr;
        colCorr = temp;
    end    
end

if ind==2
    posCorr = wc2;
    
    rowCorr = wc3-wc2; rowCorr = rowCorr ./ norm(rowCorr);
    colCorr = wc1-wc2; colCorr = colCorr ./ norm(colCorr);
    
    cost = dot(rowVectorPre, rowCorr) + dot(colVectorPre, colCorr);
    if ( cost < 0.1 )
        temp = rowCorr;
        rowCorr = colCorr;
        colCorr = temp;
    end    
end

if ind==3
    posCorr = wc3;
    
    rowCorr = wc4-wc3; rowCorr = rowCorr ./ norm(rowCorr);
    colCorr = wc2-wc3; colCorr = colCorr ./ norm(colCorr);
    
    cost = dot(rowVectorPre, rowCorr) + dot(colVectorPre, colCorr);
    if ( cost < 0.1 )
        temp = rowCorr;
        rowCorr = colCorr;
        colCorr = temp;
    end
end

if ind==4
    posCorr = wc4;
    
    rowCorr = wc1-wc4;
    colCorr = wc3-wc4;
    
    cost = dot(rowVectorPre, rowCorr) + dot(colVectorPre, colCorr);
    if ( cost < 0.1 )
        temp = rowCorr;
        rowCorr = colCorr;
        colCorr = temp;
    end    
end

rowCorr = rowCorr ./ norm(rowCorr);
colCorr = colCorr ./ norm(colCorr);

dataCorr = data;

headerCorr = CreateFtkHeaderInfo(data, pixelSpacing);
headerCorr.positionPatient = posCorr;
headerCorr.orientationPatient(1,:) = rowCorr;
headerCorr.orientationPatient(2,:) = colCorr;
headerCorr.orientationPatient(3,:) = cross(headerCorr.orientationPatient(1,:), headerCorr.orientationPatient(2,:));

for y=1:row
    for x=1:col        
        c = [x-1 y-1 0];
        [wc(1), wc(2), wc(3)] = Image2WorldMrFtk(headerCorr, x-1, y-1, 0);        
        [cPre(1), cPre(2), cPre(3)] = World2ImageMrFtk(header, wc(1), wc(2), wc(3));
        cPre = round(cPre+1);
        if ( cPre(1)>col )
            cPre(1) = col;
        end
        
        if ( cPre(2)>row )
            cPre(2) = row;
        end
        dataCorr(y, x) = data(cPre(2), cPre(1));
    end
end


    