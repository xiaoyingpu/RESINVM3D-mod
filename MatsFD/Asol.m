function[V] = Asol(MTX,u,tol)
% [V] = Asol(MTX,u,tol)
% Calculates the 3D Dc res forward modeling matrix times a vector
% it solves Av = u (in this notation) or v = A^-1*u

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2006

nit = 500;
I = speye(size(MTX.FOR));
 
A = MTX.FOR;
[n,n] = size(MTX.FOR);


fprintf('           A sol    ');
%calculate the potentials for each source term -  sources in a q = array(nx*ny*nz,# 0f source dipoles)

if size(u,1) > size(MTX.FOR,1);
    %put the long vector into a new arrangement for solving in parts
    u = reshape(u,size(MTX.FOR,1),size(u,1)/size(MTX.FOR,1));
end    
  %
%initialize the potential vector
V = zeros(n,size(u,2)); 

for i=1:size(u,2),
   qi = u(:,i);
   qi = qi(MTX.p);
 
   v = bicgstb(A,qi, MTX.PL, MTX.PU, nit, tol);
 
   %arrange as a vector
   v = v(MTX.up);
  V(:,i) =  v;

   fprintf('.');          
end;
disp(' ');
