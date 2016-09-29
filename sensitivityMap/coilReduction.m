
function [DataDst, E_0, V, dstCha] = coilReduction(Data, dstCha, dstChaThres)
% perform the coil reduction

if ( nargin < 2 )
    dstCha = floor(size(Data,3)/2);
end

if ( nargin < 3 )
    dstChaThres = 0.1;
end

[DataDst, E_0, V] = Eigen_Channel(Data, eye(size(Data,3)));
if ( dstCha  <= 0 )
    % ind = find(E_0/E_0(end)<=dstChaThres);
    ind = find(cumsum(E_0)./sum(E_0)<dstChaThres);
    if ( ~isempty(ind) )
        dstCha = size(Data, 3) - ind(end);
    else
        dstCha = size(Data, 3);
    end
%     v = E_0(1);
%     v = v * dstChaThres;
%     ind = find(E_0<=v);
%     if ( ~isempty(ind) )
%         dstCha = size(Data, 3) - ind(end);
%     else
%         dstCha = size(Data, 3);
%     end
end

if ( numel(size(Data)) == 3 )
    DataDst = DataDst(:,:,size(Data, 3)-dstCha+1:size(Data, 3));
end

if ( numel(size(Data)) == 4 )
    DataDst = DataDst(:,:,size(Data, 3)-dstCha+1:size(Data, 3),:);
end
