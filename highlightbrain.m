function [mesh_handle,mesh]=highlightbrain(pt,ROI,colr,edgeface,plotbrain,newfig,hemi)
% pt is string patient name, such as 'EC129', or just 'MNI' for standard brain
% ROI is ROIions of interest as cells of strings, 
%     for example {'superiortemporal','parsopercularis','entorhinal'}
% colr is RGB color of faces (row 1) and edges (row 2), give RGB such as [0.5 1 0.5]
% dens is 2 values for density (transparency) of surface triangles' edge
%     and face (resp.), from 0 (transparent) to 1 (opaque), example: [0.4 0.1]
% if wholebrain=1, will also bring in entire brain (or just run getbrain.m)
% if newfig==1, generates a new figure (otherwise plots in current figure)

if strcmpi(pt,'MNI'); 
    pt='cvs_avg35_inMNI152';
    disp('it really is mni!')
end
% if any(strcmp(pt,{'EC129','EC137'})); hemi='rh'; else hemi='lh'; end
if nargin<=2; colr='g'; end
if nargin<=3; edgeface=[.4 0.1]; end
% ptdir=strcat('/Users/jonkleen/Desktop/ChangLab/Data/',pt);
if newfig; figure('color','w'); end
if plotbrain; [cortex,hemi]=getbrain4(pt,plotbrain,newfig,hemi); %c_h=ctmr_gauss_plot(cortex,[],[],hemi);
% set (c_h, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.9); hold on
end
if length(hemi)==1; hemi=[hemi 'h'];end

% expand shorthand calls
isc=0;
if ischar(ROI); ROI={ROI}; isc=1; end

z=strcmp(ROI,'tt'); if any(z); ROI{1}='temporal-transverse'; end
z=strcmp(ROI,'stg'); if any(z); ROI{1}='superiortemporal'; end
z=strcmp(ROI,'mtg'); if any(z); ROI{1}='middletemporal'; end
z=strcmp(ROI,'itg'); if any(z); ROI{1}='inferiortemporal'; end
z=strcmp(ROI,'fus'); if any(z); ROI{1}='fusiform'; end
z=strcmp(ROI,'tp'); if any(z); ROI{1}='temporalpole'; end
z=strcmp(ROI,'ph'); if any(z); ROI{1}='parahippocampal'; end
z=strcmp(ROI,'hp'); if any(z); ROI{1}=[hemi 'Hippocampus']; end
z=strcmp(ROI,'am'); if any(z); ROI{1}=[hemi 'Amygdala']; end
z=strcmp(ROI,{'ent','ento','ec'}); if any(z); ROI{1}='entorhinal'; end
z=strcmp(ROI,'pt'); if any(z); ROI{1}='parstriangularis'; end
z=strcmp(ROI,{'po','pop'}); if any(z); ROI{1}='parsopercularis'; end
z=strcmp(ROI,{'por','porb'}); if any(z); ROI{1}='parsorbitalis'; end
z=strcmp(ROI,'rmf'); if any(z); ROI{1}='rostralmiddlefrontal'; end
z=strcmp(ROI,'cmf'); if any(z); ROI{1}='caudalmiddlefrontal'; end
z=strcmp(ROI,'sf'); if any(z); ROI{1}='superiorfrontal'; end
z=strcmp(ROI,'mof'); if any(z); ROI{1}='medialorbitofrontal'; end
z=strcmp(ROI,'lof'); if any(z); ROI{1}='lateralorbitofrontal'; end
z=strcmp(ROI,'fp'); if any(z); ROI{1}='frontalpole'; end
z=strcmp(ROI,'prec'); if any(z); ROI{1}='precentral'; end
z=strcmp(ROI,'postc'); if any(z); ROI{1}='postcentral'; end
z=strcmp(ROI,'precu'); if any(z); ROI{1}='precuneus'; end
z=strcmp(ROI,'parsup'); if any(z); ROI{1}='superiorparietal'; end
z=strcmp(ROI,'parinf'); if any(z); ROI{1}='inferiorparietal'; end
z=strcmp(ROI,'sm'); if any(z); ROI{1}='supramarginal'; end
z=strcmp(ROI,''); if any(z); ROI{1}=''; end
z=strcmp(ROI,'ins'); if any(z); ROI{1}='insula'; end
z=strcmp(ROI,{'cing','cingulate'}); if any(z); ROI={'rostralanteriorcingulate','caudalanteriorcingulate','posteriorcingulate','isthmuscingulate'}; end
z=strcmp(ROI,{'mfg'}); if any(z); ROI={'rostralmiddlefrontal','caudalmiddlefrontal'}; end

mesh=make_roi_mesh(pt,hemi,ROI,'test',0);
ea = edgeface(1);
fa = edgeface(2);
hold on
mesh_handle = ctmr_gauss_plot_addl(mesh, [0 0 0], 0, hemi, 0);

set (mesh_handle, 'FaceColor', colr(1,:),'EdgeColor',colr(size(colr,1),:), 'FaceAlpha', fa,'EdgeAlpha', ea);

%if strcmpi(hemi,'lh'); lightbrain('r',.8); else lightbrain('l',.8); end
cameratoolbar('setmode',''); 
