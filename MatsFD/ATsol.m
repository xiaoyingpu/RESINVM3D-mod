function[V] = ATsol(MTX,u,tol)
% [v] = ATsol(MTX,u,tol)
% Calculates the 3D DC forward modelling matrix(transposed) times a vector
%Solves the system A^t*v = u 
% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005
%
nit = 500;
A = MTX.FOR;
[n,n] = size(MTX.FOR);

fprintf('           ATsol    ');

if size(u,1) > size(MTX.FOR,1);
    %put the long vector into a new arrangement for solving in parts
    u = reshape(u,size(MTX.FOR,1),size(u,1)/size(MTX.FOR,1));
end

V = zeros(n,size(u,2));
for i=1:size(u,2),
   
   qi = u(:,i); 
%reorder for the forward operator (remeber symcrm in updateMTX)
   qi = qi(MTX.p);

   v =bicgstb(A',qi,MTX.PU',MTX.PL',nit,tol);
  v = v(MTX.up); 
   V(:,i) =  v;
 
   fprintf('.');
end;
disp(' ');