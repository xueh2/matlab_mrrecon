function plotvec(vec, names)
% Arrowline plot of a various number of 2D or 3D vectors
% SYNTAX: plotvec(<vector>, <names>)
% Plots a various number of 3D or 2D vectors.
% <vectors> has to contain a least one vector with it's components x, y (2D) or
% x, y, z (3D).
% For more than one vector the vectors had to be given in a matrix with a
% vector in each line.
% Example:
%
% <vector> = [x1 y1 ; x2 y2]
% or
% <vector> = [x1 y1 z1; x2 y2 z2; ...]
%
% <names> is optional and can contain a one-row matrix with the name of
% each vector given in <vectors>.
% The vector name is displayed inside the plot.

% AUTHOR: Hannes Eilers
% WEB: www.private-factory.de
 
% LICENCE: Creative Commons by-sa
% (http://creativecommons.org/licenses/by-sa/3.0/deed.en)

% get matrix size
sizeVec = size(vec);
n = sizeVec(1);         % rows
m = sizeVec(2);         % columns

% check matrix size
if n < 1,
    error('there had to be a least one vector')
end

if m < 2,
    error('there had to be at least two components for each vector')
end

% adding first components
u = vec( 1 : n );
% adding second components
v = vec( n + 1 : 2*n );

% adding thridn components (if exists)
if m == 3,
    w = vec( 2*n + 1 : 3*n );
end

% creating start points for vectors
x = zeros(1, n);
y = zeros(1, n);
z = zeros(1, n);

% plot vectors
if m == 3,
    quiver3( x,y,z, u,v,w, 0 )
else
    quiver( x,y, u,v, 0 )
end

% plot vector names
if nargin > 1,
    for i = 1:n,
        txt = strcat( names(i), ' (', num2str(u(i)), ';', num2str(v(i)), ';', num2str(w(i)), ')' );
        text( u(i)+0.2, v(i)+0.2, w(i)+0.2, txt )
    end
end

disp('Thank you for using plotvec. See www.private-factory.de for more.')