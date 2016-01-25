function[ MTX] = calcWTW(MTX,wt,para)
% [WTW] = calcWTW(MTX,wt)
% Calculate WTW - the model regularization matrix
% USE: grad, kron3

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005

dx = MTX.GRID.DX;
dy = MTX.GRID.DY;
dz = MTX.GRID.DZ;
%smoothing parameters for the three directions large number promotes flatness in given direction dlz big makes things vertical
alx = para.alx;
aly = para.aly;
alz = para.alz;

%Smallness parameter
als = para.als;

nx=length(dx);
ny=length(dy);
nz=length(dz);

[G,Gx,Gy,Gz] = grad(dx,dy,dz);


%%Create a weighted smallness term 
V = spdiags(mkvc(wt), 0, nx*ny*nz, nx*ny*nz); 

%Assemble the Anisotropic gradient operrator
Gs = [alx*Gx;aly*Gy;alz*Gz];

%Weights certain points more than others
Wt = spdiags(mkvc(wt),0,nx*ny*nz,nx*ny*nz);


%assemble the 3d weighting matrix
MTX.WTW = Wt' * ( Gs' * Gs + als * V) * Wt;
MTX.WTW = MTX.WTW(find(MTX.ACTIVE(:)),find(MTX.ACTIVE(:)));

