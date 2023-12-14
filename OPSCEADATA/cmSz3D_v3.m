function [cm]=cmSz3D_v3
% Adjusted for for even spread across colors
cm=[ones(64,3); getjet]; 
cm(65:76,:)=[]; 
cm(65:80,:)=[linspace(1,cm(80,1),16); linspace(1,cm(80,2),16); linspace(1,cm(80,3),16)]'; 
scm=size(cm,1); cm([1:2:round(scm*.7) round(scm*.85):2:end],:)=[];
cm=[ones(8,3); cm]; 