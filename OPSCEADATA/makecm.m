function [cm]=makecm(c_in,nrowsout)
% jon.kleen@ucsf.edu, updated 2022
% nrowsout is the how many color steps (rows of RGB) you want
% colorsin is N by 3 RGB value columns, N=number of colors
% cm is the output colormap

if ~exist('nrowsout','var'); nrowsout=64; end
if length(nrowsout)~=1; error('need a single value for nrows'); end

ntransitions=size(c_in,1)-1;
cm=[];
for i=1:ntransitions
    cm=[cm; ...
            [linspace(c_in(i,1),c_in(i+1,1),floor(nrowsout/ntransitions)); ...
             linspace(c_in(i,2),c_in(i+1,2),floor(nrowsout/ntransitions)); ...
             linspace(c_in(i,3),c_in(i+1,3),floor(nrowsout/ntransitions)); ...
             ]' ...
             
         ];
end
