% This function resize the image using k-space zero padding method to avoid
% any extra artifact caused by interpolation.
% function b = Zero_Padding_Resize(a,x,y)
% a: 2-D image or 3-D image series.
% x: size of frequence encoding direction
% y: size of phase encoding direction

function b = Zero_Padding_Resize_NoFiltering(a, x, y)

b = a;
s = size(a);
if ( x<=s(1) & y<=s(2) )
    b = a;
    return;
end

if ( length(s) == 1 )
    sx_0 = round((s(1)+1)/2); % Find the center point in k-space
    x_0 = round((x + 1)/2); % Find the center point
    t_x = x_0 - sx_0 ;
    rfft = 1/sqrt(s(1));
    rifft = sqrt(x);
    a_fft = rfft*fftshift(fft(a));
    b_fft = zeros(x);
    b_fft( 1+t_x:s(1)+t_x ) = a_fft;
    b = rifft*ifft(ifftshift(b_fft));
    return;
end

sx_0 = round((s(1)+1)/2); % Find the center point in k-space
sy_0 = round((s(2)+1)/2);

x_0 = round((x + 1)/2); % Find the center point
y_0 = round((y + 1)/2);

t_x = x_0 - sx_0 ;
t_y = y_0 - sy_0 ;

% if length(s) > 4, 'Error! Input Must Be a 2-D, 3-D or 4-D Matrix!' , return, end
if (t_x < 0) | (t_y < 0),  'Error! Output Image Size Cannot Smaller Than Input Image Size' , return, end

rfft = 1/sqrt(s(1)*s(2));
% rifft = sqrt(x*y);
rifft = (x*y)/sqrt(s(1)*s(2));

COL = s(1);
LIN = s(2);

CHA = 1;
ACQ = 1;
SLC = 1;
PAR = 1;
ECO = 1;
PHS = 1;
REP = 1;
SET = 1;

if ( numel(s) >= 3 ) CHA = s(3); end
if ( numel(s) >= 4 ) ACQ = s(4); end
if ( numel(s) >= 5 ) SLC = s(5); end
if ( numel(s) >= 6 ) PAR = s(6); end
if ( numel(s) >= 7 ) ECO = s(7); end
if ( numel(s) >= 8 ) PHS = s(8); end
if ( numel(s) >= 9 ) REP = s(9); end
if ( numel(s) >= 10 ) SET = s(10); end

b = zeros(x, y, CHA, ACQ, SLC, PAR, ECO, PHS, REP, SET);

for set=1:SET
    for rep=1:REP
        for phs=1:PHS
            for eco=1:ECO
                for par=1:PAR
                    for slc=1:SLC
                        for acq=1:ACQ
                            for cha=1:CHA                                
                                a_fft = rfft*fftshift(fft2(a(:,:,cha, acq, slc, par, eco, phs, rep, set)));
                                b_fft = zeros(x,y);
                                b_fft( 1+t_x:s(1)+t_x, 1+t_y:s(2)+t_y ) = a_fft; %imagesc(abs(b_fft))
                                b(:,:,cha, acq, slc, par, eco, phs, rep, set) = rifft*ifft2(ifftshift(b_fft));                                
                            end
                        end
                    end
                end
            end
        end
    end
end

% if length(s) == 2,
%     a_fft = rfft*fftshift(fft2(a));
%     b_fft = zeros(x,y);
%     b_fft( 1+t_x:s(1)+t_x, 1+t_y:s(2)+t_y ) = a_fft; %imagesc(abs(b_fft))
%     b = rifft*ifft2(ifftshift(b_fft));
% elseif  length(s) == 3,
%     b = zeros(x, y, s(3));
%     for i=1:s(3)
%         a2D = a(:,:,i);
%         a_fft = rfft*fftshift(fft2(a2D));
%         b_fft = zeros(x,y);
%         b_fft( (1+t_x):(s(1)+t_x), (1+t_y):(s(2)+t_y) ) = a_fft;
%         b(:,:,i) = rifft*ifft2(ifftshift(b_fft));
%     end
% elseif length(s) == 4
%     b = zeros(x, y, s(3), s(4));
%     for f=1:s(4)
%         for i=1:s(3)
%             a2D = a(:,:,i,f);
%             a_fft = rfft*fftshift(fft2(a2D));
%             b_fft = zeros(x,y);
%             b_fft( (1+t_x):(s(1)+t_x), (1+t_y):(s(2)+t_y) ) = a_fft;
%             b(:,:,i,f) = rifft*ifft2(ifftshift(b_fft));
%         end
%     end
% elseif length(s) == 5
%     b = zeros(x, y, s(3), s(4), s(5));
%     for t=1:s(5)
%         for f=1:s(4)
%             for i=1:s(3)
%                 a2D = a(:,:,i,f,t);
%                 a_fft = rfft*fftshift(fft2(a2D));
%                 b_fft = zeros(x,y);
%                 b_fft( (1+t_x):(s(1)+t_x), (1+t_y):(s(2)+t_y) ) = a_fft;
%                 b(:,:,i,f,t) = rifft*ifft2(ifftshift(b_fft));
%             end
%         end
%     end
% elseif length(s) == 6
%     b = zeros(x, y, s(3), s(4), s(5), s(6));
%     for k=1:s(6)
%         for t=1:s(5)
%             for f=1:s(4)
%                 for i=1:s(3)
%                     a2D = a(:,:,i,f,t,k);
%                     a_fft = rfft*fftshift(fft2(a2D));
%                     b_fft = zeros(x,y);
%                     b_fft( (1+t_x):(s(1)+t_x), (1+t_y):(s(2)+t_y) ) = a_fft;
%                     b(:,:,i,f,t,k) = rifft*ifft2(ifftshift(b_fft));
%                 end
%             end
%         end
%     end
% elseif length(s) == 7
%     b = zeros(x, y, s(3), s(4), s(5), s(6), s(7));
%     for p=1:s(6)
%         for k=1:s(6)
%             for t=1:s(5)
%                 for f=1:s(4)
%                     for i=1:s(3)
%                         a2D = a(:,:,i,f,t,k,p);
%                         a_fft = rfft*fftshift(fft2(a2D));
%                         b_fft = zeros(x,y);
%                         b_fft( (1+t_x):(s(1)+t_x), (1+t_y):(s(2)+t_y) ) = a_fft;
%                         b(:,:,i,f,t,k,p) = rifft*ifft2(ifftshift(b_fft));
%                     end
%                 end
%             end
%         end
%     end
% else
%     return
% end
