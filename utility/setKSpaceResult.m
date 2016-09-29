function kspace = setKSpaceResult(res, kspace, acq, slc, par, eco, rep, phs, set)
% set the kspace result

dstCHA = size(res, 3);
srcCHA = size(kspace, 3);

s = size(kspace);
if ( dstCHA ~= srcCHA )
    kspace = zeros([s(1) s(2) dstCHA s(4:end)]);
end

kspace(:,:,:,acq, slc, par, eco, rep, phs, set) = res;