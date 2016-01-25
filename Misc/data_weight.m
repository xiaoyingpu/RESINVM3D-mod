function [WdW,dtw] = data_weight(MTX,para);

%This function creates the data weighting matrix
%uses simple weighting based on the SD of the measurements

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005


%Set a threshold for the max error (percent)
maxerr = para.maxerr;
%Weights are percent standard deviation 
%get the absolute error for each datum, plus the epsilon

dtw = ((MTX.DTW/100).*abs(MTX.dobs)+para.e);
dtw = 1./dtw;

%Normalize things to 1
dtw = dtw./max(dtw);

%zero values at max error;
I = find(MTX.DTW >maxerr);
dtw(I) = 0;

WdW = spdiags(dtw, 0, length(MTX.dobs), length(MTX.dobs));
