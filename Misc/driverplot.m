function driverplot(xc,MTX);
%Plots Synthetic vs Inverse solution for the Driver files: DCdriverS and
%DCdriverBH

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky 
% Last update, August 2005.

dx = MTX.GRID.DX; nx = length(dx);
dy = MTX.GRID.DY; ny = length(dy);
dz = MTX.GRID.DZ; nz = length(dz);

[xi,yc,zc] = cell_centre(dx,dy,dz);
%%%Note xi is in the y position bc we have a left hand coord system (z down
%%%is positive)
[X,Y,Z] = meshgrid(yc,xi,zc);
figure
subplot(1,2,1);
load SynModel

slice(X,Y,Z,reshape(exp(-synmodel),nx,ny,nz),xi(ceil(nx/2)),yc(ceil(ny/2)),zc(1));
set(gca,'Zdir','reverse');
caxis([0 max(exp(-synmodel))]);
axis tight
shading flat
view([-149 5])
box on
title('Synthetic Model');
%%Reverse labeling due to LH coord system!
xlabel('Y position (m)');
ylabel('X position (m)');
zlabel('Z position (m)');


%%now make the second plot
subplot(1,2,2);


slice(X,Y,Z,reshape(exp(-xc),nx,ny,nz),xi(ceil(nx/2)),yc(ceil(ny/2)),zc(1));
set(gca,'Zdir','reverse');
caxis([0 max(exp(-synmodel))]);
axis tight
shading flat
view([-149 5])
box on
title('Inverted Model');
%%Reverse labeling due to LH coord system!
xlabel('Y position (m)');
ylabel('X position (m)');
zlabel('Z position (m)');

h= colorbar; 
a = get(h,'Ylabel');
set(a,'String','Resistivity (Ohm-m)');
