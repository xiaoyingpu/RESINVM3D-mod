function[S,Sx,Sy,Sz] = massf(s,dx,dy,dz)
% [S] = massf(s,dx,dy,dz)
%Creates a diagonal matrix that contains the harmonically averaged
%conductivity values of the faces of each cell
% s - resistivity structure dimensions (nx,ny,nz)
%dx,dy,dz are vectors containing the cell widths in the x y and z
%directions, respectively

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

% Pad for Neumann BC
%s = outint(s);

%dx = [dx(1); dx; dx(end)]; 
%dy = [dy(1); dy; dy(end)];
%dz = [dz(1); dz; dz(end)];

% Number the phi grid
np = (Nx+2)*(Ny+2)*(Nz+2); 
GRDp = reshape(1:1:np, (Nx+2),(Ny+2),(Nz+2));

% Number the Ax grid
nax = (Nx+1)*(Ny+2)*(Nz+2); 
GRDax = reshape(1:1:nax, (Nx+1),(Ny+2),(Nz+2));

% Number the Ay grid
nay = (Nx+2)*(Ny+1)*(Nz+2); 
GRDay = reshape(1:1:nay, (Nx+2),(Ny+1),(Nz+2));

% Number the Az grid
naz = (Nx+2)*(Ny+2)*(Nz+1); 
GRDaz = reshape(1:1:naz, (Nx+2),(Ny+2),(Nz+1));

% Generates the 3D grid
ex = ones(Nx+2,1);
ey = ones(Ny+2,1);
ez = ones(Nz+2,1);

Dx = kron3(dx,ey,ez);
Dy = kron3(ex,dy,ez);
Dz = kron3(ex,ey,dz);

dV = Dx.*Dy.*Dz;

%%%% Generate x coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l = 2:Nx+2; j= 1:Ny+2; k = 1:Nz+2;

% Avarage rho on x face
rhof = (dV(l,j,k).*s(l,j,k) + dV(l-1,j,k).*s(l-1,j,k))/2;

dVf = (dV(l,j,k) + dV(l-1,j,k))/2; 

rhof = rhof./dVf;
                         
lx = []; jx = []; kx = []; rx = [];

%% Coef (i,j,k)
lx = mkvc(GRDax); 
jx = mkvc(GRDax);
kx = mkvc(rhof);

Sx = sparse(lx,jx,kx);

%%%% Generate y coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l = 1:Nx+2; j= 2:Ny+2; k = 1:Nz+2;

% Avarage rho on y face
rhof = (dV(l,j-1,k).*s(l,j-1,k) + dV(l,j,k).*s(l,j,k))/2;

dVf = (dV(l,j-1,k) + dV(l,j,k))/2;

rhof = rhof./dVf;

ly = []; jy = []; ky = []; 

%% Coef (i,j,k)
ly = mkvc(GRDay);
jy = mkvc(GRDay);
ky = mkvc(rhof);

Sy = sparse(ly,jy,ky);

%%%% Generate z coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


l = 1:Nx+2; j= 1:Ny+2; k = 2:Nz+2;

% Avarage rho on z face
rhof = (dV(l,j,k-1).*s(l,j,k-1) + dV(l,j,k).*s(l,j,k))/2;

dVf = (dV(l,j,k-1) + dV(l,j,k))/2;

rhof = rhof./dVf;

lz = []; jz = []; kz = [];

%% Coef (i,j,k)
lz = mkvc(GRDaz);
jz = mkvc(GRDaz);
kz = mkvc(rhof);

Sz = sparse(lz,jz,kz);
%%%% Assemble Matrix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Oxy = sparse(nax,nay);
Oxz = sparse(nax,naz);
Oyz = sparse(nay,naz);

S = [Sx,   Oxy,   Oxz; ...
     Oxy', Sy,    Oyz; ...
     Oxz', Oyz',  Sz];
