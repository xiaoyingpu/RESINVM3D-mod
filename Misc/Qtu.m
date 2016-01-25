function[u] = Qtu(Q,d,srcnum)
% u = Qtu(Q,d)
%Calculates the solution to: u = Q'*d;

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005

u = [];
for i = 1:max(srcnum)
    j = find(srcnum == i);
    uu = Q(j,:)'*d(j);
    u = [u,uu];
end
