function[D] = div(dx,dy,dz)
% [D] = div(dx,dy,dz)
%Creates the 3D finite volume divergence operator
%operator is set up to handle variable grid discratization
%dx,dy,dz are vectors containing the cell widths in the x y and z
%directions, respectively
%
% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2006

dx = shiftdim(dx);
dy = shiftdim(dy);
dz = shiftdim(dz);

Nx = length(dx)-2;
Ny = length(dy)-2;
Nz = length(dz)-2;

%dx = [dx(1);dx;dx(end)];
%dy = [dy(1);dy;dy(end)];
%dz = [dz(1);dz;dz(end)];

% Number the phi grid 
np = (Nx+2)*(Ny+2)*(Nz+2);
GRDp = reshape(1:1:np,Nx+2,Ny+2,Nz+2);


% Number the Ax grid
nax = (Nx+1)*(Ny+2)*(Nz+2); 
GRDax = reshape(1:1:nax, (Nx+1),(Ny+2),(Nz+2));


% Number the Ay grid
nay = (Nx+2)*(Ny+1)*(Nz+2); 
GRDay = reshape(1:1:nay, (Nx+2),(Ny+1),(Nz+2));

% Number the Az grid
naz = (Nx+2)*(Ny+2)*(Nz+1); 
GRDaz = reshape(1:1:naz, (Nx+2),(Ny+2),(Nz+1));

% Generates the grid
ex = ones(Nx+2,1);
ey = ones(Ny+2,1);
ez = ones(Nz+2,1);

Dx = kron3(dx,ey,ez);
Dy = kron3(ex,dy,ez);
Dz = kron3(ex,ey,dz);



%%%%   Generate d/dx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lx = []; jx = []; kx = [];

% Entries (l,j,k)

lx = mkvc(GRDp(2:end-1,:,:));
jx = mkvc(GRDax(1:end-1,:,:));
kx = mkvc(-1./Dx(2:end-1,:,:));

% Entries (l+1,j,k)

lx = [lx;lx];
jx = [jx;mkvc(GRDax(2:end,:,:))];
kx = [kx;-kx];

% BC at x = 0

lx = [lx;mkvc(GRDp(1,:,:))];
jx = [jx;mkvc(GRDax(1,:,:))];
kx = [kx;mkvc(1./Dx(1,:,:))];

% BC at x = end

lx = [lx;mkvc(GRDp(end,:,:))];
jx = [jx;mkvc(GRDax(end,:,:))];
kx = [kx;mkvc(-1./Dx(end,:,:))];

%%%%   Generate d/dy  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ly = []; jy = []; ky = [];

% Entries (l,j,k)

ly = mkvc(GRDp(:,2:end-1,:));
jy = mkvc(GRDay(:,1:end-1,:));
ky = mkvc(-1./Dy(:,2:end-1,:));

% Entries (l+1,j,k)

ly = [ly;ly];
jy = [jy;mkvc(GRDay(:,2:end,:))];
ky = [ky;-ky];

% BC on y = 0
ly = [ly; mkvc(GRDp(:,1,:))];
jy = [jy; mkvc(GRDay(:,1,:))];
ky = [ky; mkvc(1./Dy(:,1,:))];

% BC on y = end
ly = [ly; mkvc(GRDp(:,end,:))];
jy = [jy; mkvc(GRDay(:,end,:))];
ky = [ky; mkvc(-1./Dy(:,end,:))];

%%%%   Generate d/dz  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lz = []; jz = []; kz = [];

% Entries (l,j,k)

lz = mkvc(GRDp(:,:,2:end-1));
jz = mkvc(GRDaz(:,:,1:end-1));
kz = mkvc(-1./Dz(:,:,2:end-1));

% Entries (l+1,j,k)

lz = [lz;lz];
jz = [jz;mkvc(GRDaz(:,:,2:end))];
kz = [kz;-kz];

% BC on z = 0
lz = [lz; mkvc(GRDp(:,:,1))];
jz = [jz; mkvc(GRDaz(:,:,1))];
kz = [kz; mkvc(1./Dz(:,:,1))];

% BC on z = end
lz = [lz; mkvc(GRDp(:,:,end))];
jz = [jz; mkvc(GRDaz(:,:,end))];
kz = [kz; mkvc(-1./Dz(:,:,end))];

%%%%%%%% Generate the div %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dx = sparse(lx,jx,kx,np,nax);
Dy = sparse(ly,jy,ky,np,nay);
Dz = sparse(lz,jz,kz,np,naz);

D = [Dx,Dy,Dz];
