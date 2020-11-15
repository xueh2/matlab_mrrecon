function [Fxyz] = function_create3Dmodel(MODEL)

% Function generates model for a 3D object.
% 
% SYNTAX
%   [Fxyz, N] = function_create3Dmodel(MODEL)
%
% DESCRIPTION
%   MODEL = 0; corresponds to 128x128x3 3D object which consits of 
%              'Cameraman', 'Lena', and 'Testpat' strata consiquently. 
%   MODEL = 1; corresponds to 64x64x64 3D object which consists of 
%              5 nonoverlapping spheres of different colors forming
%              'V'-shape.
%   
%   The result is Fxyz;
%
% REMARKS
%   For more details read section 9.6 'Three dimensional inverse' of the 
%   book, p. 260.
%
% Dmitriy Paliy, Tampere University of Technology, 
% Updated 31-01-2008
% dmitriy.paliy@tut.fi

if nargin<1, MODEL=1, end;

if MODEL==1,
    N = 128; % size of the object NxNxN
    
    Fxyz = zeros(N,N,N);
    
    X = [0:1/(N-1):1];
    Y = [0:1/(N-1):1];
    Z = [0:1/(N-1):1];
    
    [gX,gY,gZ] = meshgrid(X,Y,Z);

    NumberOfCells = 5;

    sX = [14:24:128]./128;
    sY = [14:24:128]./128;
    sZ = [110 62 14 62 110]./128;
    sR = [12 12 12 12 12]./128;
    sColor = [.5 .6 .7 .8 .9];
    
    for CellNumber = 1:NumberOfCells,
        Fxyz = Fxyz + sColor(CellNumber).*(sqrt( (gX - sX(CellNumber)).^2 + (gY - sY(CellNumber)).^2 + (gZ - sZ(CellNumber)).^2) <= sR(CellNumber));
    end;
elseif MODEL==2,
    forum = double(imread('forum.png'));
    forum = forum(1:2:end,1:2:end);
    image = forum;
    
    %     Cameraman = double(imread('image_Cameraman256.png'));
    %     image = Cameraman;
    
    Fxyz(:,:,1) = image;
    Fxyz = Fxyz./max(Fxyz(:));
else
    Lena = double(imread('image_Lena256.png'));
    Cameraman = double(imread('image_Cameraman256.png'));
    Testpat1 = double(imread('image_Testpat256.png'));
    %     Lena = Lena(1:2:end,1:2:end);
    %     Cameraman = Cameraman(1:2:end,1:2:end);
    %     Testpat1 = Testpat1(1:2:end,1:2:end);

    %     N = size(Lena,1);
    %     % N = 128; % size of the object NxNx3
    %
    %     Fxyz = zeros(N,N,3);
    
    Fxyz(:,:,1) = Cameraman;
    Fxyz(:,:,2) = Lena;
    Fxyz(:,:,3) = Testpat1;

    Fxyz = Fxyz./max(Fxyz(:));
end;
