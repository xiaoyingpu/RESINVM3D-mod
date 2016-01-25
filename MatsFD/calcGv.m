function[gu] = calcGv(m,u,MTX,v)
%
% [G] = calcG(m,u,MTX,v)
%
% calculate G = partial (A(m)u) / partial m
%
% calls massf, div, grad, dMassf_dm
%
% for cell center unknowns with Harmonic averaging
% assume m  = log(sigma)

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2006


dx = MTX.GRID.DX;
dy = MTX.GRID.DY;
dz = MTX.GRID.DZ;

nx = length(dx);
ny = length(dy);
nz = length(dz);

n = (nx)*(ny)*(nz);


mp = MTX.mref(:); 
mp(find(MTX.ACTIVE)) = m;
mp = reshape(mp,nx,ny,nz);

% harmonic averaged e^m into S
S = massf(exp(-mp),dx,dy,dz); 
S = spdiags(1./diag(S),0,size(S,1),size(S,2));
D = div(dx,dy,dz);
Grd = grad(dx,dy,dz);

% Obtain G. Note Gc does not depend on the model m
 Dm = spdiags(exp(-mkvc(m)),0,length(m),length(m));
% 
gu = zeros(n,size(u,2));


if size(v,2) == size(u,2);
    for i=1:size(u,2); % to number of src

        Gc = dMassf_dm(Grd*u(:,i),dx,dy,dz);
        Gc = Gc(:,find(MTX.ACTIVE));
        %calculate G*v for each source term
        gui = D*(S^2 * (Gc*(Dm*v(:,i))));;
        gu(:,i) = gui;
    end
elseif size(v,2) < size(u,2);
    for i=1:size(u,2); % to number of src
       Gc = dMassf_dm(Grd*u(:,i),dx,dy,dz);
       %calculate G*v for each source term
       Gc = Gc(:,find(MTX.ACTIVE));
       gui = D*(S^2 * (Gc*(Dm*v)));
       gu(:,i) = gui;
   end;
end;

