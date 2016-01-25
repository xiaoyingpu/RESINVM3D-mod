function [MTX,U,FOR] = generateMTX(para);
%GenerateMTX(para);
%%This function creates the data structure needed for InvMain.
%%All files must be in the .txt format *note Matlab Version 7 files cannot
%%be read by earlier versions unless saved using the -v6 flag

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky, Last update, Sept 2006

%GrdFile,SrcFile,RecFile,DataFile,MrefFile, ModeWtFile ,ActiveCells
%%Load the data files
dx =dlmread(para.dx);
dy =dlmread(para.dy);
dz =dlmread(para.dz);
srcloc =dlmread(para.SrcFile);
recloc = dlmread(para.RecFile);
data = dlmread(para.DataFile);
mref = dlmread(para.MrefFile);
model_weight = dlmread(para.ModeWtFile);
active = dlmread(para.ActiveCells);

%Create the grid related structure
MTX.GRID.DX = dx(:);
MTX.GRID.DY = dy(:);
MTX.GRID.DZ = dz(:);

%Create the RHS of the forward problem
%src_loc needs to have dimensions [length(data), 6];
%The first 3 columns are the locations of the positive electrode, columns
%4-6 are the x,y,z locations of the negative electrode. 
%if an electrode can be both positive and negative - swap locations in the
%src file and change the sign on the corresponding data point - otherwise
%you will generate extra src terms and things will run slower.

[srcnum,srcterm] = sort_src(srcloc);
[MTX.RHS] = calcRHS(MTX,srcterm);
[MTX.SRCNUM,J] = sort(srcnum);

%Sort the data to match the sources
recloc = recloc(J,1:end);
data = data(J,1:end);

%Save J for use when creating synthetic data
%%save J J;

%%Now set up the Q matrix for all reciever locations

%%first electrode
MTX.OBS = interpmat_N(dx,dy,dz,recloc(:,1),recloc(:,2),recloc(:,3)); %%%%%
%% See if it is a dipole survey - if not the other electrode is assumed to
%% be at infinity
try
MTX.OBS = MTX.OBS - interpmat_N(dx,dy,dz,recloc(:,4),recloc(:,5),recloc(:,6));
end

%%Create the data related components;
%the observed data
MTX.dobs = data(1:end,1);
%Data weighting - percent standard deviation of data
MTX.DTW = data(1:end,2);

%Model related components;
%reference model
MTX.mref = mref;
%list active cells, 1 for active; 0 for inactive
MTX.ACTIVE = active;
%Model weighting vector, for penalizing certain components of the model
MTX.wt = model_weight;
%Create the regularization matrix
[MTX] = calcWTW(MTX,model_weight,para);

if strcmp((para.BC), 'yes');
    [MTX,U,FOR] = boundary_correction(MTX,para,srcterm);
end



