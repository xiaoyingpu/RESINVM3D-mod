%Driver file for a demo inversion of surface resistivity data;
%Copyright (c) 2007 by the Society of Exploration Geophysicists.
%For more information, go to http://software.seg.org/2007/0001 .
%You must read and accept usage terms at:
%http://software.seg.org/disclaimer.txt before use.
% 
%Revision history:
%Original SEG version by Adam Pidlisecky, Oct 2006.


% Adding the path for the inversion code
addpath LinSolve Optimize Misc MatsFD ModelData ModelData/Surface


para.BETA = 5; %Reg parameter use para.BETA =[]; for auto selection
para.maxit = 5; %number of GN iterations
para.tol = 1e-3;  %Tolerance of gradient reduction (stop point)

% Inner iteration pars
para.intol = 2e-4;     %  tol for inexact newton solver (ipcg)
para.inintol = 1e-9;     %  tol for the forward and adjoint problems 
para.ininintol = 1e-6;   %  tol for the inner solution in the ipcg
para.init = 3;           %  number of ipcg iterations
para.ilutol = 1e-4;     %ilu preconditioner tolerance, reduce to 1e-3 if you run into memory issues
para.alp=1e-4;          % Parameter for line search

%Smoothing parameters larger alx,aly,alz promotes more smoothing in the
%corresponding direction; als is the smallness factor
para.alx = 1; %x flatness
para.aly = 1;  %y flatness
para.alz = 1;  %z flatness
para.als = 1e-2;

%%Data Weighting Params
%minimum epsilion for data weighting so small values don't get too much weight;
para.e = 0.1; %%must be greater than smallest abs(dobs) 
para.maxerr = 10; %max Standard Dev that is acceptable

% Input files - comma delimited ASCII

para.dx = 'dx.txt'; %must contain 1 vectors labeled that designate cell dimensions
para.dy = 'dy.txt'; %must contain 1 vectors labeled that designate cell dimensions
para.dz = 'dz.txt'; %must contain 1 vectors labeled that designate cell dimensions
para.SrcFile = 'SrcFile.txt'; %a [ndata, 6] matrix containing the x,y,z coords of the +ve and -ve source locations
para.RecFile ='RecFile.txt'; %for a pole-dipole survey this is an [ndata, 3] matrix
                         %containing the x,y,z coords of the +ve receiver
                         %for a dipole-dipole survey this is an [ndata, 6] matrix
                         %containing the x,y,z coords of the +ve and -ve receiver locations
para.DataFile = 'DataFile.txt';%a [ndata, 2] matrix containing the observed data in 1st column and standard deviation (percent)
                           % in the second column  Matrix must be labeled data;  
para.MrefFile = 'MrefFile.txt'; %a vector, length [1, number_model] containing the reference model conductivity structure
                           
para.ModeWtFile = 'ModeWtFile.txt'; %a vector, length [1, number_model] containing the relative model weights
                                
para.ActiveCells = 'ActiveCells.txt'; %a vector, length [1, number_model] containing 0 for inactive cells and 1 for active cells
                             

%FLag for applying Boundary Condition/Sigularity correction 'no' is default
para.BC = 'no'; %'yes' 'no'

 try
%     %See if the file exists
     load MTXS.mat;
     disp('Loading MTX structure');
 catch
    %if not create the file
    disp('Creating MTX structure');
    MTXS = generateMTX(para);
    save MTXS.mat MTXS;
 end;

MTXS, para
[xc,d,MTXS,para] = InvMainN(MTXS,para);

driverplot(xc,MTXS);

rmpath LinSolve Optimize Misc MatsFD ModelData ModelData/Surface