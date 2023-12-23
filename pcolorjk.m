function [m]=pcolorjk(m)
% adds a row and column of nans so that pcolor can be used w/o losing row
% and column from visualization


m=[m  nan(size(m,1),1)];
m=[m; nan(1,size(m,2))];


pcolor(m); 


