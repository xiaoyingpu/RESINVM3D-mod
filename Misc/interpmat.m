function [Q] = interpmat(dx,dy,dz,xr,yr,zr)
% [Q] = interpmat(dx,dy,dz,xr,yr,zr)
% Interpolation matrices for potential field 
%xr,yr,zr are vectors containing the locations of the sources or receivers
%in real space. the code assume coords are centered at 0,0 in the xy plane,
%and z is 0 at the surface and positive down (LH system)
%
%Will blow up if receiver is located on the outside cells
%
%
% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005


dx = shiftdim(dx);
dy = shiftdim(dy);
dz = shiftdim(dz);

%build the 3d grid - numbered from 0 to maximum extent
z(1) = 0; for i=1:length(dz); z(i+1) = z(i)+dz(i); end;
x(1) = 0; for i=1:length(dx); x(i+1) = x(i)+dx(i); end;
y(1) = 0; for i=1:length(dy); y(i+1) = y(i)+dy(i); end;

%Center the grid about zero
x = shiftdim(x) - max(x)/2;
y = shiftdim(y) - max(y)/2;

%z = shiftdim(z) - max(z)/2;
%Set surface to Z = 0
z= shiftdim(z);

%find the cell centers
xc = x(1:end-1) +dx/2;
yc = y(1:end-1) +dy/2;
zc = z(1:end-1) +dz/2;

%call linear interp scheme below
[Q] = linint(xc,yc,zc,xr,yr,zr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Q] = linint(x,y,z,xr,yr,zr)
%
% This function does local linear interpolation
% computed for each receiver point in turn
%
% calls mkvc
%
% [Q] = linint(x,y,z,xr,yr,zr)
% Interpolation matrix 
%

 
nx = length(x) ;
ny = length(y) ;
nz = length(z) ;

np = length(xr);

Q = sparse(np,nx*ny*nz);

for i = 1:np,

       % fprintf('Point %d\n',i); 

        [dd,im] = min(abs(xr(i)-x));
        if  xr(i) - x(im) >= 0,  % Point on the left 
                 ind_x(1) = im;
                 ind_x(2) = im+1;
        elseif  xr(i) - x(im) < 0,  % Point on the right
                 ind_x(1) = im-1; 
                 ind_x(2) = im;
       end;
       dx(1) = xr(i) - x(ind_x(1));
       dx(2) = x(ind_x(2)) - xr(i);

        [dd,im] = min(abs(yr(i) - y)) ; 
       if  yr(i) - y(im) >= 0,     % Point on the left
                 ind_y(1) = im;
                 ind_y(2) = im+1;
       elseif  yr(i) -y(im) < 0,  % Point on the right
                 ind_y(1) = im-1;
                 ind_y(2) = im;
       end;
       dy(1) = yr(i) - y(ind_y(1));
       dy(2) = y(ind_y(2)) - yr(i);

        [dd,im] = min(abs(zr(i) - z));
        if  zr(i) -z(im) >= 0,  % Point on the left
                 ind_z(1) = im;
                 ind_z(2) = im+1;
       elseif  zr(i) -z(im) < 0,  % Point on the right
                 ind_z(1) = im-1;
                 ind_z(2) = im;
       end;
       dz(1) = zr(i) - z(ind_z(1)); 
       dz(2) = z(ind_z(2)) - zr(i);      

       dv = (x(ind_x(2)) - x(ind_x(1))) * (y(ind_y(2)) - y(ind_y(1))) * ...
            (z(ind_z(2)) - z(ind_z(1)));

      Dx =  (x(ind_x(2)) - x(ind_x(1)));
      Dy =  (y(ind_y(2)) - y(ind_y(1)));
      Dz =  (z(ind_z(2)) - z(ind_z(1)));  

      
      % Get the row in the matrix
      v = zeros(nx, ny,nz);

      v( ind_x(1),  ind_y(1),  ind_z(1)) = (1-dx(1)/Dx)*(1-dy(1)/Dy)*(1-dz(1)/Dz);
      v( ind_x(1),  ind_y(2),  ind_z(1)) = (1-dx(1)/Dx)*(1-dy(2)/Dy)*(1-dz(1)/Dz);
      v( ind_x(2),  ind_y(1),  ind_z(1)) = (1-dx(2)/Dx)*(1-dy(1)/Dy)*(1-dz(1)/Dz);
      v( ind_x(2),  ind_y(2),  ind_z(1)) = (1-dx(2)/Dx)*(1-dy(2)/Dy)*(1-dz(1)/Dz);
      v( ind_x(1),  ind_y(1),  ind_z(2)) = (1-dx(1)/Dx)*(1-dy(1)/Dy)*(1-dz(2)/Dz);
      v( ind_x(1),  ind_y(2),  ind_z(2)) = (1-dx(1)/Dx)*(1-dy(2)/Dy)*(1-dz(2)/Dz);
      v( ind_x(2),  ind_y(1),  ind_z(2)) = (1-dx(2)/Dx)*(1-dy(1)/Dy)*(1-dz(2)/Dz);
      v( ind_x(2),  ind_y(2),  ind_z(2)) = (1-dx(2)/Dx)*(1-dy(2)/Dy)*(1-dz(2)/Dz);

     
      Q(i,:) = mkvc(v)';

end;


