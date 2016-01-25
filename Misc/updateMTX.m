function[MTX] = updateMTX(MTX,m,ilutol)
% [MTX] = updateMTX(MTX,m,ilutol)
% Update the MTX structure for a new model
% m is the new model

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005

dx = MTX.GRID.DX; nx = length(dx);
dy = MTX.GRID.DY; ny = length(dy);
dz = MTX.GRID.DZ; nz = length(dz);
m = m(:);

mp = MTX.mref(:); 
mp(find(MTX.ACTIVE)) = m;


MTX.mc = m;
mp = reshape(mp,nx,ny,nz);
%Generate the conductivity structure

Msig = massf(exp(-mp),dx,dy,dz);
Msig = spdiags(1./diag(Msig),0,size(Msig,1),size(Msig,2));

% %generate the operators
% MTX.Div = div(dx,dy,dz);
% MTX.Grad = grad(dx,dy,dz);
% MTX.MSIG = Msig;

%make the forward modeling matrix
MTX.FOR = div(dx,dy,dz) * Msig * grad(dx,dy,dz);

%This pins the corner of the potential field by a constant
 B = sparse(nx*ny*nz,nx*ny*nz); 
 B(1,1) = 1/dx(1)/dy(1)/dz(1);

 %Adjust MTX.FOR
MTX.FOR(1,:) = MTX.FOR(1,:) + B(1,:);

%move things closer to the diag by permutation (symrcm)
%for faster linear solve
p = symrcm(MTX.FOR);
[ii,up] = sort(p);
MTX.FOR = MTX.FOR(p,p);

%calculate the preconditioner
% TODO: luinc deprecated 
% [PL,PU] = luinc(MTX.FOR, ilutol);
setup.droptol = ilutol;
[PL,PU] = ilu(MTX.FOR, setup);
MTX.PL = PL;
MTX.PU = PU;
MTX.up = up;
MTX.p = p;