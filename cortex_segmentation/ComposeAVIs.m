
function ComposeAVIs(avifilenames, header, outfileName, outfileAVIName)

% if ( isempty(dir(outfileAVIName))==0 )
%     disp(['Already exsiting: ' outfileAVIName]);
%     return;
% end

if ( iscell(avifilenames)==0 )
    return;
end

if ( numel(avifilenames)==0 )
    return;
end

N = numel(avifilenames);

m = aviread(avifilenames{1});

num = size(m, 2);

colormap = m(1).colormap;
[H, W, D] = size(m(1).cdata);

% for k=1:N
%     M = CheckMovieSize2(M, H, W);
% end
    

% different situation
colnum = 3;
rownum = ceil(N/3);

if ( D==1 )
    M = zeros(rownum*H, colnum*W, num, 'uint8')+128;

    for kk=1:N
        
        rowHead = H *floor(kk/(colnum+1)) + 1
        colHead = W *mod(kk-1, colnum)+ 1    
        
        avifilenames{kk}
        
        if ( isempty(dir(avifilenames{kk}) ) )
            continue;
        end
        m = aviread(avifilenames{kk});       
        m = CheckMovieSize2(m, H, W);
        
        for tt=1:num
        
            M( rowHead: rowHead+H-1, colHead:colHead+W-1 , tt) = m(tt).cdata;    
        
        end
    
    end

%     for tt=1:num
%         
%         M2(tt) = im2frame(uint8(M(:, :, tt)), gray(256));
%     
%     end    

    header.xsize = colnum*W;
    header.ysize = rownum*H;

    data2avi(M, outfileAVIName, 1);
%     SaveAnalyze(uint32(M), header, outfileName, 'Grey');

    clear M
%     hdr2avi(outfileName, outfileAVIName, 1);
end

if ( D==3 )
    M = zeros(rownum*H, colnum*W, D, num, 'uint8');

    for kk=1:N
        avifilenames{kk}
        if ( isempty(dir(avifilenames{kk}) ) )
            continue;
        end
        m = aviread(avifilenames{kk});
        m = CheckMovieSize2(m, H, W);
        
        rowHead = H *floor(kk/(colnum+1)) + 1
        colHead = W *mod(kk-1, colnum)+ 1    
    
        for tt=1:num
        
            M( rowHead: rowHead+H-1, colHead:colHead+W-1 , 1:3, tt) = double(m(tt).cdata);    
        
        end
    
    end

    for tt=1:num
        
        M2(tt) = im2frame(M(:, :, 1:3, tt));
    
    end    

    FramesRate = 8;
    global frameRate
    if ( frameRate ~= -1 )
        FramesRate = frameRate;
    end

    movie2avi(M2, outfileAVIName, 'Compression', 'None', 'fps', FramesRate, 'quality', 100);
end