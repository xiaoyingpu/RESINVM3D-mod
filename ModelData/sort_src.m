function [srcnum,srcterm] = sort_src(src_loc);
%[srcnum,srcterm] = sort_src(src_loc);
%%This function sorts a list of source terms and creates an index term for
%%use in creating the MTX structure

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky, Last update, July 2006

%initialize the sort index
[sort_temp, sort_ind] = sortrows(src_loc,[1:size(src_loc,2)]);

%%Now find out how many independent current sources there are
sort_l = sort_temp(1:end-1,:)-sort_temp(2:end,:);

I = find(sum(abs(sort_l),2));
%now we make the source term;

srcterm = sort_temp([I; I(end)+1 ],:);

%%Now create the source numbering sequence;
srcnum = zeros(length(src_loc),1);


for i = 1:length(I);
    if i == 1;
        srcnum(1:I(1)) = i;
    else
       srcnum(I(i-1)+1:I(i)) = i;
    end;
end;

srcnum(I(end)+1:end) = length(I)+1;


%%Reshuffle the srcnum term so that it matches the data.

srcnum(sort_ind) = srcnum;;
