function[S,Sx,Sy,Sz] = dMassf_dm(u,dx,dy,dz)
% [S,Sx,Sy,Sz] = dMassf_dm(u,dx,dy,dz)
% Differentiate the face matrix
% creates the a matrix containing the derivative of the face matrix (massf)
% u is a 3-D matrix, dimensions(length(dx),length(dy),length(dz));
%dx,dy,dz are vectors containing the cell widths in the x y and z
%directions, respectively
% Calls kron3, mkvc

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


% Number the phi grid
np = (Nx+2)*(Ny+2)*(Nz+2); 
GRDs = reshape(1:1:np, (Nx+2),(Ny+2),(Nz+2));

% Number the Ax grid
nax = (Nx+1)*(Ny+2)*(Nz+2); 
GRDx = reshape(1:1:nax, (Nx+1),(Ny+2),(Nz+2));

% Number the Ay grid
nay = (Nx+2)*(Ny+1)*(Nz+2); 
GRDy = reshape(1:1:nay, (Nx+2),(Ny+1),(Nz+2));

% Number the Az grid
naz = (Nx+2)*(Ny+2)*(Nz+1); 
GRDz = reshape(1:1:naz, (Nx+2),(Ny+2),(Nz+1));

% Generates the 3D grid
ex = ones(Nx+2,1);
ey = ones(Ny+2,1);
ez = ones(Nz+2,1);

Dx = kron3(dx,ey,ez);
Dy = kron3(ex,dy,ez);
Dz = kron3(ex,ey,dz);

dV = Dx.*Dy.*Dz;

ux = u(1:nax);
uy = u(nax+1:nax+nay);
uz = u(nax+nay+1:nax+nay+naz);

Ux = reshape(ux,Nx+1,Ny+2,Nz+2);
Uy = reshape(uy,Nx+2,Ny+1,Nz+2);
Uz = reshape(uz,Nx+2,Ny+2,Nz+1);

%%%% Generate x coefficienticients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = 2:Nx+2; j = 1:Ny+2; k = 1:Nz+2;

dVf = (dV(l,j,k) + dV(l-1,j,k))/2;

C1 = GRDx*0; C2 = C1;

C1 = dV(l-1,j,k)./dVf/2;
C2 = dV(l,j,k)./dVf/2;
C1 = C1.*Ux;
C2 = C2.*Ux;

lx = mkvc(GRDx);
jx = mkvc(GRDs(l-1,j,k)); 
kx = mkvc(C1);

lx = [lx;lx];
jx = [jx; mkvc(GRDs(l,j,k))];
kx = [kx;mkvc(C2)];

inz = find(jx);
lx = lx(inz); jx = jx(inz); kx = kx(inz);

Sx = sparse(lx,jx,kx,nax,np);

%%%% Generate y coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = 1:Nx+2; j= 2:Ny+2; k = 1:Nz+2;

dVf = (dV(l,j,k) + dV(l,j-1,k))/2;

C1 = GRDy*0; C2 = C1;

C1 = dV(l,j-1,k)./dVf/2;
C2 = dV(l,j,k)./dVf/2;
C1 = C1.*Uy;
C2 = C2.*Uy;

ly = mkvc(GRDy);
jy = mkvc(GRDs(l,j-1,k)); 
ky = mkvc(C1);

ly = [ly;ly];
jy = [jy; mkvc(GRDs(l,j,k))];
ky = [ky;mkvc(C2)];

inz = find(jy);
ly = ly(inz); jy = jy(inz); ky = ky(inz);

Sy = sparse(ly,jy,ky,nay,np);

%%%% Generate z coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = 1:Nx+2; j= 1:Ny+2; k = 2:Nz+2;

dVf = (dV(l,j,k) + dV(l,j,k-1))/2;

C1 = GRDz*0; C2 = C1;

C1 = dV(l,j,k-1)./dVf/2;
C2 = dV(l,j,k)./dVf/2;
C1 = C1.*Uz;
C2 = C2.*Uz;

lz = mkvc(GRDz);
jz = mkvc(GRDs(l,j,k-1)); 
kz = mkvc(C1);

lz = [lz;lz];
jz = [jz; mkvc(GRDs(l,j,k))];
kz = [kz;mkvc(C2)];

inz = find(jz);
lz = lz(inz); jz = jz(inz); kz = kz(inz);

Sz = sparse(lz,jz,kz,naz,np);

%%%% Assemble Matrix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Oxy = sparse(nax,nay);
Oxz = sparse(nax,naz);
Oyz = sparse(nay,naz);

S = [Sx;Sy;Sz];








