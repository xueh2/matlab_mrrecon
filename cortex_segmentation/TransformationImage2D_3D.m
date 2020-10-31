
function sourceT = TransformationImage2D_3D(target, source, header, dx, dy)

sourceT = zeros(size(source), 'uint32');

for k=1:header.zsize
    k
    for j=1:header.ysize
        for i=1:header.xsize
            ti = dx(j, i, k)+i;
            tj = dy(j, i, k)+j;
            value = interp2(source(:,:,k), ti, tj, 'linear');
            sourceT(j, i, k) = value;
        end
    end
end