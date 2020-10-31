
function ComposeAVIs_MN_GUI(M_Frames, outfileAVIName, N_col, FrameRate)
% M rows and N cols

N = numel(M_Frames);

H = -1;
W = -1;
SliceNum = -1;
for kk=1:N
    
    [h, w, D] = size(M_Frames{kk}(1).cdata);
    slicenum = size(M_Frames{kk}, 2);
    
    if ( h>H )
        H = h;
    end
    
    if ( w>W )
        W = w;
    end
    
    if ( slicenum>SliceNum )
        SliceNum = slicenum;
    end
end

if ( mod(N, N_col)==0 )
    M_row = N / N_col;
else
    M_row = 1 + floor(N / N_col);
end

MI = zeros(M_row*H, N_col*W, 3, SliceNum, 'uint8');

for kk=1:N
    
    m = M_Frames{kk};
    slicenum = size(m, 2);
    [h, w, D] = size(m(1).cdata);
    
    rowHead = H *floor(kk/(N_col+1)) + 1
    colHead = W *mod(kk-1, N_col)+ 1    

    for tt=1:slicenum
        cdata = zeros( H, W, 3 , 'uint8');
        
        if ( D==3 )
            cdata(1:h, 1:w, 1:3) = m(tt).cdata;
        end
        
        if ( D==1 )
            cdata(1:h, 1:w, 1) = m(tt).cdata;
            cdata(1:h, 1:w, 2) = m(tt).cdata;
            cdata(1:h, 1:w, 3) = m(tt).cdata;
        end
        
        MI( rowHead: rowHead+H-1, colHead:colHead+W-1 , 1:3, tt) = cdata;    

    end

end

for tt=1:SliceNum

    M2(tt) = im2frame(MI(:, :, 1:3, tt));

end    

movie2avi(M2, outfileAVIName, 'Compression', 'None', 'fps', FrameRate, 'quality', 100);
