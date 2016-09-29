function dstRecon = PFConvRecon2D(src, dst, rangeFE, rangePE, kSize, kernel, rawKernel)
% src, dst: src and dst kspace
% rangeFE, rangePE: sampled region
% kernel: kSize(1)*kSize(2)*srcCHA*dstCHA
% dstRecon: reconed dst kspace

[sx,sy,srcCHA, N] = size(src);
[sx,sy,dstCHA, N] = size(dst);

% rawKernel = reshape(kernel, [kSize(1)*kSize(2)*srcCHA dstCHA]);

D = [];
R = [];

% R = D*rawkernel
halfKx = floor(kSize(1)/2);
halfKy = floor(kSize(2)/2);

% composite the src data matrix

tmp = zeros(N, kSize(1)*kSize(2)*srcCHA);

indRow = 1;
for pe=1:sy
    for fe=1:sx
        
        if ( (fe<=rangeFE(2)) & (fe>=rangeFE(1)) & (pe<=rangePE(2)) & (pe>=rangePE(1)) )
            continue;
        end
        
        indCol = 1;
        for s=1:srcCHA
            for y=-halfKy:halfKy
                for x=-halfKx:halfKx

                    px = fe + x;
                    py = pe + y;

                    if ( px < 1 )
                        px = px + sx;
                    end

                    if ( px > sx )
                        px = px - sx;
                    end

                    if ( py < 1 )
                        py = py + sy;
                    end

                    if ( py > sy )
                        py = py - sy;
                    end
                    
                    tmp(:, indCol) = reshape(src(px, py, s, :), [N 1]);

                    indCol = indCol + 1;
                end
            end
        end
        
        D = [D; tmp];             
        indRow = indRow + N;
    end
end

R = D * rawKernel;

dstRecon = dst;

indRow = 1;
for pe=1:sy
    for fe=1:sx
        
        if ( (fe<=rangeFE(2)) & (fe>=rangeFE(1)) & (pe<=rangePE(2)) & (pe>=rangePE(1)) )
            continue;
        end
        
        for n=1:N
            dstRecon(fe, pe, :, n) = R(indRow+n-1, :);
        end
        
        indRow = indRow + N;
    end
end

