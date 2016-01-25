function[v] = mkvc(A)
% v = mkvc(A)
%reshapes a 3D matrix A, into a vector.
%note this function can be replace with the call A = A(:)

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Eldad Haber
% Last update, July 2001

v = reshape(A,size(A,1)*size(A,2)*size(A,3),1);