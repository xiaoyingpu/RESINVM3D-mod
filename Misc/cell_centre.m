function [xc,yc,zc] = cell_centre(dx,dy,dz);
%%[xc,yc,zc] = cell_centre(dx,dy,dz);
%%Finds the LH coordsystem cartesian coords of the cell centers;

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005

%calculate the cartesian cell centered grid for inversion
dx = mkvc(dx);
dy = mkvc(dy);
dz = mkvc(dz);
%build the 3d grid - numbered from 0 to maximum extent
z = [0; cumsum(dz)];
x = [0; cumsum(dx)];
y = [0; cumsum(dy)];

%Center the grid about zero
x = shiftdim(x) - max(x)/2;
y = shiftdim(y) - max(y)/2;

%z = shiftdim(z) - max(z)/2;
%Set surface to Z = 0
z= shiftdim(z);

%find the cell centers
xc = x(1:end-1) +dx/2;
yc = y(1:end-1) +dy/2;
zc = z(1:end-1) +dz/2;