function[q] = calcRHS(MTX,src);
% [q] = calcRHS(MTX,src)
% Generate the 3D DC resistivity forward modeling RHS
%%Uses MTX from generate MTX
%%src is an n*6 matrix with the locations of the positive and negative
%%source locations as follows [xpos ypos zpos xneg yneg zneg]; n is the
%%number of RHS terms

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky, Last update, August 2006


dx = MTX.GRID.DX;
dy = MTX.GRID.DY;
dz = MTX.GRID.DZ;

nx = length(dx);
ny = length(dy);
nz = length(dz);
%%generate a vector that contains each cell volume
dv = mkvc(kron3(dx,dy,dz));
%%%%%%%% Loop over sources - where sources are a dipole %%%%%%%%%%%
% allocate space
nh = nx*ny*nz;

fprintf('   Calc RHS (for source)   ');

%if there is more than one source   
for k=1:size(src,1);
    %interpolate the location of the sources to the nearest cell nodes
    Q = interpmat_N(dx,dy,dz,src(k,1),src(k,2),src(k,3)); %%%%%
    Q = Q - interpmat_N(dx,dy,dz,src(k,4),src(k,5),src(k,6));
    
    q(:,k) = mkvc(Q)./dv; 
    %fprintf('.')   
end;
disp('Done ');