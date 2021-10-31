
function [starting_points, ending_points] = CCMS_Contour(I, level, connectivity, displayflag, contourColor, contourWidth)
% [starting_points, ending_points] = CCMS_Contour(I, level, connectivity, displayflag, contourColor, contourWidth)
% according to the connectivity consistency to determine the contours for a 2D data
% starting_points, ending_points store the starting and ending coordinates for every short lines
% the coordinates: x from left to right (column), y from top to bottom (row), leftup point is (1,1)

[ysize, xsize] = size(I);
starting_points = zeros(xsize*ysize, 2);
ending_points =  zeros(xsize*ysize, 2);
num = 1;

index = find(I==level);
if ( isempty(index) == 0 )
    I(index) = level-1e-4;
end

ps = [0 0];
pe = [0 0];
for j=1:ysize-1
    for i=1:xsize-1
        % for points
        p1 = [i j];
        v1 = I(j, i);
        L1 = v1>level;
        
        p2 = [i+1 j];
        v2 = I(j, i+1);
        L2 = v2>level;
        
        p3 = [i+1 j+1];
        v3 = I(j+1, i+1);
        L3 = v3>level;
        
        p4 = [i j+1];
        v4 = I(j+1, i);
        L4 = v4>level;
        
        sumL = L1+L2+L3+L4;
        
        if (sumL==0)
            % outside
            continue;
        end
        
        if (sumL==4)
            % inside
            continue;
        end
        
        if ((sumL==1) | (sumL==3))
            % no ambigity, add one lines
            
            if ( sumL==3 )
                L1=~L1;
                L2=~L2;
                L3=~L3;
                L4=~L4;
            end
            
            if ( L1 )
                % p1-p4, p1-p2
                ps = p1_p4(p1,v1,p4,v4,level);
                pe = p1_p2(p1,v1,p2,v2,level);

                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue;
            end
            
            if ( L2 )
                % p1-p2, p3-p2
                ps = p1_p2(p1,v1,p2,v2,level);
                pe = p3_p2(p3,v3,p2,v2,level);

                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue;
            end

            if ( L3 )
                % p3-p4, p3-p2
                ps = p3_p4(p3,v3,p4,v4,level);
                pe = p3_p2(p3,v3,p2,v2,level);

                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue;
            end

            if ( L4 )
                % p1-p4, p3-p4
                ps = p1_p4(p1,v1,p4,v4,level);
                pe = p3_p4(p3,v3,p4,v4,level);

                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue;
            end
        end
                
        if (sumL==2)
            
            if ( L1&L2 )
                
                % p1-p4, p3-p2
                ps = p1_p4(p1,v1,p4,v4,level);
                pe = p3_p2(p3,v3,p2,v2,level);
                
                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue; 
            end

            if ( L1&L4 )
                
                % p1-p2, p3-p4
                ps = p1_p2(p1,v1,p2,v2,level);
                pe = p3_p4(p3,v3,p4,v4,level);
                
                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue;
            end
            
            if ( L2&L3 )
                % p1-p2, p3-p4
                ps = p1_p2(p1,v1,p2,v2,level);
                pe = p3_p4(p3,v3,p4,v4,level);
                
                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue;
            end
            
            if ( L3&L4 )
                % p1-p4, p3-p2
                ps = p1_p4(p1,v1,p4,v4,level);
                pe = p3_p2(p3,v3,p2,v2,level);

                starting_points(num,:) = ps;
                ending_points(num,:) = pe;
                num = num+1;
                continue;
            end           
            
            if ( L2&L4 )
                
                if ( connectivity==4 )
                    % p4-p1, p4-p3
                    ps = p1_p4(p1,v1,p4,v4,level);
                    pe = p3_p4(p3,v3,p4,v4,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;
                    
                    % p2-p1, p2-p3
                    ps = p1_p2(p1,v1,p2,v2,level);
                    pe = p3_p2(p3,v3,p2,v2,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;
                    continue;
                end
                
                if ( connectivity==8 )
                    % p4-p1, p2-p1
                    ps = p1_p4(p1,v1,p4,v4,level);
                    pe = p1_p2(p1,v1,p2,v2,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;

                    % p4-p3, p2-p3
                    ps = p3_p4(p3,v3,p4,v4,level);
                    pe = p3_p2(p3,v3,p2,v2,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;
                    continue;
                end
                
            end
            
            if ( L1&L3 )
                                
                if ( connectivity==4 )
                    % p1-p4, p1-p2
                    ps = p1_p4(p1,v1,p4,v4,level);
                    pe = p1_p2(p1,v1,p2,v2,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;

                    % p3-p4, p3-p2
                    ps = p3_p4(p3,v3,p4,v4,level);
                    pe = p3_p2(p3,v3,p2,v2,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;
                    continue;
                end
                
                if ( connectivity==8 )
                    % p1-p4, p3-p4
                    ps = p1_p4(p1,v1,p4,v4,level);
                    pe = p3_p4(p3,v3,p4,v4,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;
                    
                    % p1-p2, p3-p2
                    ps = p1_p2(p1,v1,p2,v2,level);
                    pe = p3_p2(p3,v3,p2,v2,level);
                    
                    starting_points(num,:) = ps;
                    ending_points(num,:) = pe;
                    num = num+1;
                    continue;
                end
            end
            
        end
    end
end

starting_points = starting_points(1:num-1,:);
ending_points = ending_points(1:num-1,:);

if nargin < 6, contourWidth = 2; end
if nargin < 5, contourColor = [1 0 0]; end
if nargin < 4, displayflag = 0; end

if (displayflag)
    f=figure;
    hold on
    imshow(I,[]);
    for tt=1:num-1
       line([starting_points(tt,1) ending_points(tt,1)],[starting_points(tt,2) ending_points(tt,2)], ...
           'LineWidth', contourWidth, 'Color', contourColor); 
    end
    hold off
end


%-------------------------------------------------------------
function p = p1_p4(p1,v1,p4,v4,level)
i = p1(1);
j = p1(2);
y = LinearInterpolation(j+1, v4, j, v1, level);
p = [i y];

%-------------------------------------------------------------
function p = p1_p2(p1,v1,p2,v2,level)
i = p1(1);
j = p1(2);
x = LinearInterpolation(i, v1, i+1, v2, level);
p = [x j];

%-------------------------------------------------------------
%-------------------------------------------------------------

function p = p3_p2(p3,v3,p2,v2,level)
i = p3(1)-1;
j = p3(2)-1;
y = LinearInterpolation(j+1, v3, j, v2, level);
p = [i+1 y];

%-------------------------------------------------------------
function p = p3_p4(p3,v3,p4,v4,level)
i = p3(1)-1;
j = p3(2)-1;
x = LinearInterpolation(i+1, v3, i, v4, level);
p = [x j+1];

%-------------------------------------------------------------
%-------------------------------------------------------------

function x = LinearInterpolation(x1, v1, x2, v2, v)
x = (x2*(v-v1)+x1*(v2-v))/(v2-v1);