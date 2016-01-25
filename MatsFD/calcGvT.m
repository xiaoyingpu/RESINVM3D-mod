function[gu] = calcGvT(m,u,MTX,v)
% [G] = calcG(m,u,nx,ny,nz)
%
% calculate G = partial (A(m)u) / partial m
%
% calls massf, div, grad, dMassf_dm
%
% for cell center unknowns with Harmonic averaging
% assume sigma = exp(m)

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
% G = D*S^2 * Gc*Dm;
 Dm = spdiags(exp(-mkvc(m(:))),0,length(m(:)),length(m(:)));
        

gu = zeros(length(m),size(u,2));
%gu=[];
for i=1:size(u,2); % to number of src

    Gc = dMassf_dm(Grd*u(:,i),dx,dy,dz);
    Gc = Gc(:,find(MTX.ACTIVE(:)));
    
    %calculate G*v for each source term

    gui =  Dm'*(Gc'*((S^2)'*(D'*v(:,i))));
    gu(:,i) =gui;
    
end
gu = sum(gu,2);
