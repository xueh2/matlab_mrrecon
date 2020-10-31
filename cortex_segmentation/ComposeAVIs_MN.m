
function ComposeAVIs_MN(avifilenames, header, outfileName, outfileAVIName, N_col)
% M rows and N cols

if ( isempty(dir(outfileAVIName))==0 )
    disp(['Already exsiting: ' outfileAVIName]);
    return;
end

N = numel(avifilenames);

for i=1:N
    p = dir(avifilenames{i});
    if ( isempty(p) )
        disp(['NOT exsiting: ' avifilenames{i}]);
        return;
    end
end

m = aviread(avifilenames{1});

num = size(m, 2);

colormap = m(1).colormap;
[H, W, D] = size(m(1).cdata);

if ( mod(N, N_col)==0 )
    M_row = N / N_col;
else
    M_row = 1 + floor(N / N_col);
end

if ( D==1 )
    MI = zeros(M_row*H, N_col*W, num)+128;

    for kk=1:N
        
        rowHead = H *floor(kk/(N_col+1)) + 1
        colHead = W *mod(kk-1, N_col)+ 1    
        
        avifilenames{kk}
        
        if ( isempty(dir(avifilenames{kk}) ) )
            continue;
        end
        m = aviread(avifilenames{kk});       
    
        for tt=1:num
        
            MI( rowHead: rowHead+H-1, colHead:colHead+W-1 , tt) = m(tt).cdata;    
        
        end
    
    end

    for tt=1:num
        
        M2(tt) = im2frame(uint8(MI(:, :, tt)), gray(256));
    
    end    

    header.xsize = 3*W;
    header.ysize = 2*H;

    SaveAnalyze(uint32(MI), header, outfileName, 'Grey');

    hdr2avi(outfileName, outfileAVIName, 1);
end

if ( D==3 )
    MI = zeros(M_row*H, N_col*W, D, num);

    for kk=1:N
        avifilenames{kk}
        if ( isempty(dir(avifilenames{kk}) ) )
            continue;
        end
        m = aviread(avifilenames{kk});
    
        rowHead = H *floor(kk/(N_col+1)) + 1
        colHead = W *mod(kk-1, N_col)+ 1    
    
        for tt=1:num
        
            MI( rowHead: rowHead+H-1, colHead:colHead+W-1 , 1:3, tt) = double(m(tt).cdata);    
        
        end
    
    end

    for tt=1:num
        
        M2(tt) = im2frame(MI(:, :, 1:3, tt)/255);
    
    end    

    movie2avi(M2, outfileAVIName, 'Compression', 'None', 'fps', 8, 'quality', 100);
end