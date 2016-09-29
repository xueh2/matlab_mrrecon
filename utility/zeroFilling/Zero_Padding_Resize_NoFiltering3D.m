% This function resize the image using k-space zero padding method to avoid
% any extra artifact caused by interpolation.
% function b = Zero_Padding_Resize3D(a,x,y,z)
% a: 5-D image [COL LIN ACQ SLC PAR]
% x: size of frequence encoding direction
% y: size of phase encoding direction
% z: size of partition direction

function b = Zero_Padding_Resize_NoFiltering3D(a, x, y, z)

b = a;
s = size(a);
if ( x==s(1) & y==s(2) & z==s(5) )
    b = a;
    return;
end

sx_0 = round((s(1)+1)/2); % Find the center point in k-space
sy_0 = round((s(2)+1)/2);
sz_0 = round((s(5)+1)/2);

x_0 = round((x + 1)/2); % Find the center point
y_0 = round((y + 1)/2);
z_0 = round((z + 1)/2);

t_x = x_0 - sx_0 ;
t_y = y_0 - sy_0 ;
t_z = z_0 - sz_0 ;

% if length(s) > 4, 'Error! Input Must Be a 2-D, 3-D or 4-D Matrix!' , return, end
if (t_x < 0) | (t_y < 0) | (t_z < 0)
    return;
end

rfft = 1/sqrt(s(1)*s(2)*s(5));
rifft = (x*y*z)/sqrt(s(1)*s(2)*s(5));

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

CHA = 1;
if ( numel(s) >= 3 ) ACQ = s(3); end
if ( numel(s) >= 4 ) SLC = s(4); end
if ( numel(s) >= 5 ) PAR = s(5); end
if ( numel(s) >= 6 ) ECO = s(6); end
if ( numel(s) >= 7 ) PHS = s(7); end
if ( numel(s) >= 8 ) REP = s(8); end
if ( numel(s) >= 9 ) SET = s(9); end

b = zeros(x, y, ACQ, SLC, z, ECO, PHS, REP, SET);

for set=1:SET
    for rep=1:REP
        for phs=1:PHS
            for eco=1:ECO
                for slc=1:SLC
                    for acq=1:ACQ
                        %for cha=1:CHA                                                            
                            a_fft = rfft*fftshift(fftn(squeeze(a(:,:,acq, slc, :, eco, phs, rep, set))));
                            b_fft = zeros(x,y,z);
                            b_fft( 1+t_x:s(1)+t_x, 1+t_y:s(2)+t_y, 1+t_z:s(5)+t_z ) = a_fft;
                            b(:,:,acq, slc, :, eco, phs, rep, set) = rifft*ifftn(ifftshift(b_fft));                                
                        %end
                    end
                end
            end
        end
    end
end
