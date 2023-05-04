% function max_avg_MNI(sz_nns_mat,sz_w8s_mat,szxyz_mat,mni_xyz_cell, npt,hem)


mni_mesh_path = '/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/OPSCEADATA/avg_change_folders/cvs_avg35_inMNI152/Imaging/Meshes/';

cd('/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/OPSCEADATA/avg_change_folders/cvs_avg35_inMNI152/Imaging/Meshes/');


npt=length(manual_ptsz);

%4/28 TO DO: 
%   1.  INDEX E- WITHIN 5MM OF EACH VERTEX INTO VECTOR (USE "FIND" TO INDEX)
%   2.  DO LOGICAL INDEX ON THE MNI XYZ COORDINATES (
%   3.  FIND CORRESPONDING SEIZURE XYZ COORDINATE/E- W8S 
%   4.  GET MAX VALUE OF W8S OF E- WITHIN VECTOR




%%


%INPUT HERE
dst_radius = 10; % minimum distance in mm from electrode to each vertex
hem='l';


%START
load([mni_mesh_path 'cvs_avg35_inMNI152_' hem 'h_pial.mat']); %load right mni vertices
mni_h = cortex.vert;
nvert=size(mni_h,1);

max_w8s_rad=nan(nvert,npt);
numelecclose=nan(nvert,npt);
tic
for sz_i = 1:npt % loop through each patient seizure
        em = mni_xyz_cell{sz_i}; % warped MNI electrodes from specific patient
        w8s = sz_w8s_mat{1,sz_i}(sz_nns_mat{1,sz_i}); % patient-specific electrode vector of ll mean difference before and after seizure symptom
    parfor v_m = 1:nvert % loop through each vertex in mni h pial

        if size(em,1)~=length(w8s); error('em does not match w8s'); end
        dst=nan(size(w8s));
           for e = 1:size(em,1) % loop through each electrode
              dst(e,1) = distance3D(em(e,:),mni_h(v_m,:)); % distance in mm from vertex (v_m) to each warped mni electrode coordinate (e) 
           end
              elecclose=dst<dst_radius;
              if any(elecclose)
                  numelecclose(v_m,sz_i)=sum(elecclose); % number of electrodes near that vertex for that patient
                  max_w8s_rad(v_m,sz_i)=max(w8s(elecclose)); % max LL change value across the electrodes near that vertex for that patient
              end
    end
end

haselecclose=numelecclose>0; %logical index of whether a patient has at least one electrode near that vertex





%%
figure('color','w','position',[230 171 1440 796])
minnumpts=2;
nptAtVerts=sum(haselecclose,2);
nptPositive=sum(max_w8s_rad>0,2); 
percentPositive=nptPositive./nptAtVerts; 

% plot3(mni_h(nptAtVerts<minnumpts,1),mni_h(nptAtVerts<minnumpts,2),mni_h(nptAtVerts<minnumpts,3),'.','color',[0 0 0],'markersize',3);

getbrain4_ns('MNI','',1,0,hem)
shading interp
axis equal; axis off; set(gca,'clipping','off'); lightsout; litebrain(hem,1)
hold on
    pause(1)
for v_m = 1:nvert
    if nptAtVerts(v_m)>=minnumpts
        plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'.','color',[1 -percentPositive(v_m)+1 -percentPositive(v_m)+1],'markersize',15);
        disp([num2str(round(v_m/nvert*100,2)) '%'])
    end
end

%%
disp('you did it!')


return

% This is all in MNI coordinates (MNI mesh for vertices, and the MNI warped electrodes, which we can call clinical_elecs_all_warped.mat file)
%   
% Goal is to end up with a vector of values called "VectorIndex" that is the same length as the number of vertices on the MNI brain. (there are thousands!)
% 
% 
% 1. For each vertex, treat it like an entire structure (like pre central gyrus):
%     a. DONE v_idx: For each patient, find electrodes (MNI coordinates) within     5mm of that vertex (MNI) 
%     b. DONE w8s_rad_cell: get corresponding weights of those electrodes
%     c. DONE max_w8s_rad_cell: take the max value of those electrodes (again, for each patient). 
%        For example, if 6 patients have electrodes within 5m of that vertex, you would end up with 6 values, regardless of whether the patients have 1, 2 or 8 electrodes, etc.
%     d. DONE: mean_v: Record the mean values of those values from "a" above as the value for that entry of VectorIndex.
%     end
% end

% add ll meandiff weights  (instead of max distance get max weights) -- need to import!
% DONE: make elecbook a cell array before importing it 
% DONE: change radius to 5mm instead of 10 
% par for loop of mni vertices as main loop, pt loop embedded within (look up how to do par for loop)
    % NOTE: MUST HAVE PT AS OUTER LOOP IN ORDER TO INDEX THROUGH CONVERTED SZ --> MNI VERTEX CELL ARRAY TO DETERMINE LATERALITY
% DONE:need to make sure ON THE CORRECT SIDE OF BRAIN 
% proportion? 

% INITIALIZE VARIABLES
% 
% load([mni_mesh_path 'cvs_avg35_inMNI152_' hem 'h_pial.mat']); %load right mni vertices
% mni_rh = cortex.vert;
% Rcortex = cortex;
% clear cortex
% 
% load([mni_mesh_path 'cvs_avg35_inMNI152_' hem 'h_pial.mat']); %load left mni vertices
% mni_lh = cortex.vert;
% Lcortex = cortex;
% clear cortex
% 
% cd('/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA')
% 
% mni_max = max(length(mni_rh),length(mni_lh)); % max number of MNI vertices of right and left hem
% dst_cell = cell(mni_max,npt); %cell array to collect distance in mm of each vertex to each electrode 
% logic_rad_cell = cell(mni_max,npt); %cell array to collect if vertex has electrode 5mm away from it (1) or not (0)
% 
% v_idx = cell(mni_max,npt); %cell array to collect index of all electrodes 5mm away from each vertex
% 
% w8s_rad_cell = cell(mni_max,npt); %cell array to collect weight of each electrode on each MNI vertex
% 
% max_w8s_rad_cell_r = cell(mni_max,npt); %cell array to collect  max weight across electrodes on each right MNI vertex
% max_w8s_rad_cell_l = cell(mni_max,npt); %cell array to collect  max weight across electrodes on each left MNI vertex
% 
% mean_v_r = zeros(length(mni_max),1); %vectors to collect mean of max weight across electrodes on each right MNI vertex
% mean_v_l = zeros(length(mni_max),1); %vectors to collect mean of max weight across electrodes on each left MNI vertex
% 
%  %%
% for v_m = 1:mni_max % loop through each vertex in mni h pial
%     for sz_i = 1:npt % loop through each patient seizure
% 
%         em = mni_xyz_cell{sz_i}; % warped MNI electrodes from specific patient
%         isR = nansum(em(:,1))>0; % laterality right
%         isL=isR~=1; %laterality left
% 
%         if (isR) && (v_m <= length(mni_rh))
%            for e = 1:size(em,1) % loop through each electrode
%               dst_cell{v_m,sz_i}(e,1) = distance3D(em(e,:),mni_rh(v_m,:)); % distance in mm from vertex (v_m) to each warped mni electrode coordinate (e) 
%               logic_rad_cell{v_m,sz_i}(e,1) = dst_cell{v_m,sz_i}(e,1) < dst_radius; % logical cell array of whether distance is less than dst_radius (1) or not (0)
%            end
%         elseif (isL) && (v_m <= length(mni_lh))
%            for e = 1:size(em,1) % loop through each electrode
%                dst_cell{v_m,sz_i}(e,1) = distance3D(em(e,:),mni_lh(v_m,:)); % distance in mm from vertex (v_m) to each warped mni coordinate (e) 
%                logic_rad_cell{v_m,sz_i}(e,1) = dst_cell{v_m,sz_i}(e,1) < dst_radius; % logical cell array of whether distance is less than dst_radius (1) or not (0)
%            end
%         end
%     end
% end
% 
% 
% 
% %%
% 
% 
% 
% clear v_m sz_i
% 
% f1 = figure;
% f2 = figure;
% %%
% for v_m = 1:mni_max 
%     for sz_i = 1:npt
%           w8s = sz_w8s_mat{1,sz_i}(sz_nns_mat{1,sz_i}); % patient-specific electrode vector of ll mean difference before and after seizure symptom
% 
%           em = mni_xyz_cell{sz_i}; % warped MNI electrodes from specific patient
%           isR = nansum(em(:,1))>0; % laterality right
%           isL=isR~=1; %laterality left
% 
%           if find(logic_rad_cell{v_m,sz_i})
%               v_idx{v_m,sz_i} = find(logic_rad_cell{v_m,sz_i}); %indices of vertices 5mm away from electrode
%               w8s_rad_cell{v_m,sz_i} = w8s(v_idx{v_m,sz_i}); %w8s(v_idx{v_m,s_i}) is vector of 1 weight per electrode (eg. [-.0432, .64]), max takes makes value per partient per vertex
%               if isR 
%                   max_w8s_rad_cell_r{v_m,sz_i} = max(w8s(v_idx{v_m,sz_i})); %max weight of vertex for each patient 
%               elseif isL 
%                   max_w8s_rad_cell_l{v_m,sz_i} = max(w8s(v_idx{v_m,sz_i})); %max weight of vertex for each patient 
%               end
%           end
%     end        
%     %mean_v_r(v_m,1) = mean([max_w8s_rad_cell_r{v_m,:}]); %mean vector of max weights across all patients
%     %mean_v_l(v_m,1) = mean([max_w8s_rad_cell_l{v_m,:}]); %mean vector of max weights across all patients
%     sum_v_r(v_m,1) = sum([max_w8s_rad_cell_r{v_m,:}]>0); %mean vector of max weights across all patients
%     sum_v_l(v_m,1) = sum([max_w8s_rad_cell_l{v_m,:}]>0); %mean vector of max weights across all patients
%     
%     
%     
% 
%     if v_m <= length(mni_rh)
%         figure(f1)
%         plot3(mni_rh(v_m,1),mni_rh(v_m,2),mni_rh(v_m,3),'.','color',[0 0 0]);
%         hold on
%     end
%     if v_m <= length(mni_lh)
%         figure(f2)
%         plot3(mni_lh(v_m,1),mni_lh(v_m,2),mni_lh(v_m,3),'.','color',[0 0 0]);
%         hold on
%     end
% end
% 
% axis equal
% 
% 
% cd('/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA')
% 
% % 
% %  