function [MTX,U,FOR] = boundary_correction(MTX,para,srcterm);
%%Removes singulrity and decreases boundary effects by modifying source
%%term based on the analytical solution for a homogeneous halfspace

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky, Last update, July 2006


dx = MTX.GRID.DX;
dy = MTX.GRID.DY;
dz = MTX.GRID.DZ;

nx = length(MTX.GRID.DX);
ny = length(MTX.GRID.DY);
nz = length(MTX.GRID.DZ);
%%this function 

%%First Create a grid in realspace

%build the 3d grid - numbered fromn 0 to maximum extent
z(1) = 0; for i=1:length(dz); z(i+1) = z(i)+dz(i); end;
x(1) = 0; for i=1:length(dx); x(i+1) = x(i)+dx(i); end;
y(1) = 0; for i=1:length(dy); y(i+1) = y(i)+dy(i); end;

%Center the grid about zero
x = shiftdim(x) - max(x)/2;
y = shiftdim(y) - max(y)/2;

%Set surface to Z = 0
z= shiftdim(z);

%find the cell centers
xc = x(1:end-1) +dx/2;
yc = y(1:end-1) +dy/2;
zc = z(1:end-1) +dz/2;

[X,Y,Z] = ndgrid(xc,yc,zc);

U = zeros(size(MTX.RHS));
%solve for u on this grid using average mref;

avg_cond = geomean(exp(mkvc(MTX.mref)));

%turn off warning b/c we get a divide by zero that we will fix later
warning('off');
%loop over all sources
for i = 1:size(MTX.RHS,2);
    
pve1 = ((X - srcterm(i,1)).^2 + (Y - srcterm(i,2)).^2 + (Z - srcterm(i,3)).^2).^0.5;
%norm of negative current electrode and 1st potential electrode
nve1 = ((X - srcterm(i,4)).^2 + (Y - srcterm(i,5)).^2 + (Z - srcterm(i,6)).^2).^0.5;
%norm of imaginary positive current electrode and 1st potential electrode
pveimag1 = ((X - srcterm(i,1)).^2 + (Y - srcterm(i,2)).^2 + (Z + srcterm(i,3)).^2).^0.5;
%norm of imaginary negative current electrode and 1st potential electrode
nveimag1 = ((X - srcterm(i,4)).^2 + (Y - srcterm(i,5)).^2 + (Z + srcterm(i,6)).^2).^0.5;
U(:,i) = reshape(1/(avg_cond*4*pi)*(1./pve1-1./nve1+1./pveimag1-1./nveimag1),length(mkvc(MTX.mref)),1);
end
warning('on');
%%now check for singularities due to the source being on a node
for i = 1:size(MTX.RHS,2);
    I  = find(isinf(U(:,i)));
    
    if max(size(I)) > 0;
        for j = 1:length(I);
            [a,b,c] = ind2sub([nx,ny,nz], I(j));
            %%Check to see if this a surface electrode
            if c ==1
                         U(I(j),i) = mean(U(sub2ind([nx,ny,nz],a+1,b,c),i) +  U(sub2ind([nx,ny,nz],a,b+1,c),i)...
                + U(sub2ind([nx,ny,nz],a,b,c+1),i) +  U(sub2ind([nx,ny,nz],a-1,b,c),i)...
                +  U(sub2ind([nx,ny,nz],a,b-1,c),i));
            else
            U(I(j),i) = mean(U(sub2ind([nx,ny,nz],a+1,b,c),i) +  U(sub2ind([nx,ny,nz],a,b+1,c),i)...
                + U(sub2ind([nx,ny,nz],a,b,c+1),i) +  U(sub2ind([nx,ny,nz],a-1,b,c),i)...
                +  U(sub2ind([nx,ny,nz],a,b-1,c),i) +  U(sub2ind([nx,ny,nz],a,b,c-1),i));
            end;
        end;
    end;
end;

Msig = massf(ones(nx,ny,nz)*(1./avg_cond),dx,dy,dz);
Msig = spdiags(1./diag(Msig),0,size(Msig,1),size(Msig,2));
D = div(dx,dy,dz);
G = grad(dx,dy,dz);


FOR = D * Msig* G;
 B = sparse(nx*ny*nz,nx*ny*nz); 
 B(1,1) = 1/dx(1)/dy(1)/dz(1);
FOR(1,:) = FOR(1,:) + B(1,:);
 
for i = 1:size(MTX.RHS,2);

    MTX.RHS(:,i) = FOR*U(:,i);
end;