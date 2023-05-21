function [cortexout,hem]=getbrain4_ns(pt_name,sz_name,plotbrain,newfig,hem)
% pt is patient such as EC137
% plotbrain is 1 to plot the brain or 0 to not plot it (e.g. in case you are just running getbrain.m to get the outputs)
% newfig is 1 to create a new figure for this plot
% hem is l, r, or b (left, right, bilateral)
% updated by Natalia Sucher 1/17/2023

mainpath='/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA/';      %path for OPSCEA folders



if ~exist(mainpath,'dir'); mainpath='/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA/'; end
%     if ~exist(mainpath,'dir'); mainpath='/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA'; end
    
if ~exist('plotbrain','var') || isempty(plotbrain)
    plotbrain=0;
end

if ~exist('newfig','var') || isempty(newfig)
    newfig=0; 
end

if ~exist('hem','var') || isempty(hem)
    hem='b';
end %default: bilateral

if any(strcmp(hem,{'b','lr'}))
    bilat = 1;
else
    bilat = 0;
end 

if strcmpi(pt_name,'MNI')
    pt_name='cvs_avg35_inMNI152'; 
end 

ptdir=[mainpath 'OPSCEADATA/avg_change_folders/' pt_name '/Imaging/Meshes/' pt_name];
% if ~[exist([ptdir '_lh_pial.mat'])==2]; ptdir=[mainpath(1:end-7) 'AN/DATA/' pt '/Imaging/Meshes/' pt]; end
% if ~[exist([ptdir '_lh_pial.mat'])==2]; ptdir=['/Users/jkleen/Desktop/ChangLab/Data/' pt '/Imaging/Meshes/' pt]; end
if ~[exist([ptdir '_lh_pial.mat'])==2]
    ptdir=['/nataliasucher/Desktop/UCSF/Coding/OPSCEA/OPSCEADATA/avg_change_folders/' pt_name '/Imaging/Meshes/' pt_name]; 
end

% ptdir=[mainpath 'OPSCEADATA/' pt '/Imaging/Meshes/' pt];
% if ~[exist([ptdir '_lh_pial.mat'])==2]; ptdir=[mainpath(1:end-7) 'AN/DATA/' pt '/Imaging/Meshes/' pt]; end



cortexout={};
if nargin>=2
    if plotbrain
        if newfig
            figure('color','w','Name',[pt_name '-' sz_name])
        end

        if bilat
            hold on; 
        end
        
        if bilat||strcmp(hem(1),'l') %plot left brain
            load([ptdir '_lh_pial.mat']); 
            ctmr_gauss_plot_addl(cortex,[],[],'lh'); 
            cortexout=[cortexout cortex];
%             litebrain('l',.25)
        end
            
        if bilat||strcmp(hem(1),'r') % plot right brain
            load([ptdir '_rh_pial.mat']); 
            ctmr_gauss_plot_addl(cortex,[],[],'rh'); 
            cortexout=[cortexout cortex];
%             litebrain('r',.25)
        end

        if bilat
            litebrain('a',.25); 
        end
    elseif nargout>0
        if bilat||strcmp(hem(1),'l') %get left brain
            load([ptdir '_lh_pial.mat']); cortexout=[cortexout cortex];
        end
        if bilat||strcmp(hem(1),'r') % get right brain
            load([ptdir '_rh_pial.mat']); cortexout=[cortexout cortex];
        end
    end
end

l = 0;

