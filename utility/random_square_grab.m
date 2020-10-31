function [grab_image,direction,x_loc,y_loc] = random_square_grab(InputImage, SizeGrabImage)
% Returns a cropped image of the input image. The cropped image is
% generated from a random pixel location and with a random orientation.

% There are 8 orientations of the cropped image as per the 8-pixel
% neighborhood of any pixel

%Input:
% InputImage=The matrix representing the color/grayscale image
% SizeGrabImage= Dimensions of the cropped grayscale image

% Output:
% grab_image= Cropped gray-scale image
% direction= Orientation of the output image (1,2,..,8)
% x_loc, y_loc= X and Y co-ordinate of the pixel representating the pixel
% (1,1) of the cropped image

% Usage: (To crop a randomly orientated 100 X 100 image)
% I=imread('cameraman.tif');
% [grab_image,direction,x_loc,y_loc] = random_square_grab(I,100);

% Contact Author:
% Soumyabrata Dev
% http://www3.ntu.edu.sg/home2012/soumyabr001/

if (length(size(InputImage))==3)    % If it is a RGB Image, convert it to a RGB image
    I=rgb2gray(InputImage);
    [Y,X]=size(I);
else
    I=InputImage;
    [Y,X]=size(I);
end

grab_image=zeros(SizeGrabImage,SizeGrabImage); % Size of the grab image


while (1)
    
    try
    
    direction=round(1+7*rand);  % Orientation of the grabed image. There are 8 orientations as per 8-neighborhood pixels in an image

    switch direction
        case 1
            INV=I'; % Rows and columns swapped in order to make the image matrix in similar manner as the co-ordinate system
            x_start=round((X-SizeGrabImage)*rand);   % X-pixel location of the grabbed image
            y_start=round((Y-SizeGrabImage)*rand);   % Y-pixel location of the grabbed image
            for i=1:SizeGrabImage
                x=x_start;  y=y_start+i-1;
                for j=1:SizeGrabImage
                    grab_image(i,j)=INV(x,y);
                    x=x+1;
                end
            end
         
        case 2
            INV=I';
            x_start=round((X-SizeGrabImage*2)*rand);
            y_start=round(SizeGrabImage+(Y-SizeGrabImage -SizeGrabImage)*rand);
            for i=1:SizeGrabImage
                x=x_start-1+i;  y=y_start-1+i;
                    for j=1:SizeGrabImage
                        grab_image(i,j)=INV(x,y);
                        x=x+1;
                        y=y-1;
                    end
            end
    
        case 3
            INV=I';
            x_start=round((X-SizeGrabImage)*rand);  
            y_start=SizeGrabImage+round((Y-SizeGrabImage)*rand);
            for i=1:SizeGrabImage
                x=x_start+i-1;  y=y_start;
                for j=1:SizeGrabImage
                    grab_image(i,j)=INV(x,y);
                    y=y-1;
                end
            end  
        
        case 4
            INV=I';
            x_start=round(SizeGrabImage+(X-SizeGrabImage -SizeGrabImage)*rand);
            y_start=round(SizeGrabImage*2+(Y-SizeGrabImage*2)*rand);
            for i=1:SizeGrabImage
                x=x_start+i-1;  y=y_start-i+1;
                for j=1:SizeGrabImage
                    grab_image(i,j)=INV(x,y);
                    x=x-1;
                    y=y-1;
                end
            end       
        
        case 5
            INV=I';
            x_start=round(SizeGrabImage+(X-SizeGrabImage)*rand);
            y_start=round(SizeGrabImage+(Y-SizeGrabImage)*rand);
            for i=1:SizeGrabImage
                x=x_start;  y=y_start-i+1;
                for j=1:SizeGrabImage
                    grab_image(i,j)=INV(x,y);
                    x=x-1;
                end
            end        
        
        
        case 6
            INV=I';
            x_start=round(SizeGrabImage*2+(X-SizeGrabImage*2)*rand);
            y_start=round(SizeGrabImage*2+round((Y-SizeGrabImage*2)*rand));
            for i=1:SizeGrabImage
                x=x_start;  y=y_start-i+1;
                for j=1:SizeGrabImage
                    grab_image(i,j)=INV(x,y);
                    x=x-1;
                end
            end         
        
        case 7
            INV=I';
            x_start=SizeGrabImage+round((X-SizeGrabImage)*rand);  
            y_start=round((Y-SizeGrabImage)*rand);
            for i=1:SizeGrabImage
                x=x_start-i+1;  y=y_start;
                for j=1:SizeGrabImage
                    grab_image(i,j)=INV(x,y);
                    y=y+1;
                end
            end          
        
        
        case 8
            INV=I';
            x_start=round(SizeGrabImage+(X-SizeGrabImage -SizeGrabImage)*rand);
            y_start= round((Y-2*SizeGrabImage)*rand);
            for i=1:SizeGrabImage
                x=x_start-i+1;  y=y_start+i-1;
                for j=1:SizeGrabImage
                    grab_image(i,j)=INV(x,y);
                    x=x+1;
                    y=y+1;
                end
            end          
        
        otherwise
            disp ('No other direction available');
        
    end

    catch
     
        continue;
        
    end

    break;
end

%figure(1); imshow(uint8(grab_image));   % Smaller grabbed image
%figure(2); imshow(I); % Original input image

x_loc=x_start;  y_loc=y_start;