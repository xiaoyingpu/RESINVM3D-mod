function[v] = Qu(Q,u,srcnum)
% [v] = Qu(Q,u)
%projects the data onto the Q grid
%solves the linear system v = Q*u

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005

v = [];

for i = 1:size(u,2)
    %find q cells related to the source config
    j = find(srcnum == i)  ;
    vv = Q(j,:)*u(:,i);
    v = [v;vv];
end

