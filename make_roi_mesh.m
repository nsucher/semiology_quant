function mesh = make_roi_mesh(subj,hem,label_list,roi_name,debug_plot,sav)
% function mesh = make_roi_mesh(subj, hem, label_list)
%
% Input:
%   subj [str]: 'EC63', e.g.
%   hem [str]: 'lh' or 'rh'
%   label_list [cell array]: list of the labels to include in your ROI,
%                            e.g. label_list = {'precentral','postcentral'}
%   roi_name  [str]: what you want to call the ROI, e.g. 'vSMC', for file
%                   name
%   debug_plot [1 or 0]: whether to plot the ROI after it's created
%   (default  =0)
%
%   save [1 or 0]: whether to save the ROI after it's created
%   (default  =0)
%
% Output:
%   mesh structure including only vertices from the labels in label_list
%
% Format of label file is as follows:
% [vertex #]  [x coord]  [y coord]  [z coord]  [dummy value - freesurfer says "don't know" what this is]
%
if nargin <5
    debug_plot = 0;
end

fsdir =  '/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA/OPSCEADATA'; %'/Volumes/Emma/OPSCEADATA'; %'/Volumes/Sz3D/Sz3DDATA'; % Directory with your pial meshes and the labels
if exist(fsdir)==0; fsdir='/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA/OPSCEADATA'; end
subjdir = sprintf('%s/%s', fsdir, subj);
%if ~exist([subjdir , '/Imaging/Meshes']); fsdir =  '/Users/jkleen/Desktop/ChangLab/Data'; subjdir = sprintf('%s/%s', fsdir, subj); end %may be in ChangLab folder
recondir = sprintf('%s/Imaging/Meshes', subjdir); % for cvs_avg35_inMNI152

labeldir = sprintf('%s/Imaging/label/gyri', subjdir);
% outfile = sprintf('%s/%s_%s_%s_pial.mat', recondir, subj, hem, roi_name);
outfile = sprintf('%s/%s_%s_pial.mat', recondir, subj, hem, roi_name);

% Load the cortical surface mesh
% load(sprintf('%s/%s_%s_pial.mat', recondir, subj, hem));
load(sprintf('%s/%s_%s%s_pial.mat', recondir, subj, hem(1),hem(2)));


pial_surf = struct();
pial_surf.cortex = cortex;
%pial_surf = load(sprintf('%s/%s_%s_pial.mat', recondir, subj, hem));
mesh = struct();
mesh.tri = [];
mesh.vert = [];
vertnums = [];
for lab = 1:length(label_list)
    this_label=sprintf('%s/%s.%s.label', labeldir, hem, label_list{lab});
    if exist(this_label,'file')
      fid=fopen(this_label);
      C = textscan(fid, '%f %f %f %f %f','Headerlines',2);
      verts = C{1}+1;
      
      % Find the vertices and add them to the new mesh vertices
      vertnums = [verts; vertnums];
      fclose(fid);
    end
end

% Sort the vertices so they're drawn in the correct order
vertnums = sort(vertnums);
mesh.vert = pial_surf.cortex.vert(vertnums,:);

vnum_new = 1:length(vertnums); % Index of the vertex (new, relative to ofc lobe)
tri_list = [];

% Find the triangles for these vertex numbers
tri_row_inds = find(sum(ismember(pial_surf.cortex.tri, vertnums),2)==3);
tri_list = pial_surf.cortex.tri(tri_row_inds,:);

[i,j]=ismember(tri_list,vertnums);

mesh.tri=j;

if exist('sav','var') && sav;
fprintf(1,'Saving %s recon in %s\n', roi_name, outfile);
save(outfile, 'mesh');
end

if debug_plot
    ctmr_gauss_plot(mesh, [0 0 0], 0, hem);
end
