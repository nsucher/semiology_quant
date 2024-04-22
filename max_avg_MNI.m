function max_avg_MNI(sx_input, sz_nns_mat,sz_w8s_mat,mni_xyz_cell,npt,hem,dst_radius,minnumpts,opscea_path,data_path,mp_count)


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


mni_mesh_path = [data_path 'cvs_avg35_inMNI152/Imaging/Meshes/'];

cd(mni_mesh_path);


%START
load([mni_mesh_path 'cvs_avg35_inMNI152_' hem 'h_pial.mat']); %load right mni vertices
mni_h = cortex.vert;
nvert=size(mni_h,1);

cd(opscea_path)

max_w8s_rad=nan(nvert,npt);
numelecclose=nan(nvert,npt);
tic
for sz_i = 1:npt % loop through each patient seizure
        em = mni_xyz_cell{sz_i}; % warped MNI electrodes from specific patient
        if ~isempty(sz_nns_mat{1,sz_i})
            w8s = sz_w8s_mat{1,sz_i}(sz_nns_mat{1,sz_i}); % patient-specific electrode vector of ll mean difference before and after seizure symptom
        else
            w8s = []; % patient-specific electrode vector of ll mean difference before and after seizure symptom
        end
    parfor v_m = 1:nvert % loop through each vertex in mni h pial

        if size(em,1)==length(w8s)
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
end

haselecclose=numelecclose>0; %logical index of whether a patient has at least one electrode near that vertex

nptAtVerts=sum(haselecclose,2);
nptPositive=sum(max_w8s_rad>0,2); 
% 
% pink_lemonade = colormap(flipud(colormap(spring)));
% pink_lemonade = colormap(spring);

cm_percent = cbrewer2('Reds',150,'cubic'); %color map 
cm_percent = cm_percent(20:120,:);
% 
percentPositive= round(nptPositive./nptAtVerts * length(cm_percent)); 

% percentPositive= round(nptPositive./nptAtVerts * length(pink_lemonade)); 


fig_name_vert_pos = [sx_input{1,1}, ' ', num2str(minnumpts), ' patients significant positive activity change', hem];

% % PINK LEMONADE BRAIN OF % OF PATIENTS WITH POSITIVE ACTIVITY
% % 



figure('color','w','position',[230 171 1440 796],'Name', fig_name_vert_pos);
getbrain4_ns('MNI','',1,0,hem,opscea_path,data_path) %plot empty brain
shading interp
axis equal; axis off; set(gca,'clipping','off'); lightsout; litebrain(hem,.75)
hold on

pause(1)


for v_m = 1:nvert %for this specific vertex
    if nptAtVerts(v_m)>=minnumpts %if there are at least minnumpts number of patients with electrodes dst_radius away from the vertex
        if percentPositive(v_m) == 0
            percentPositive(v_m) = 1;
        end
%         plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'.','color',pink_lemonade(percentPositive(v_m),:),'markersize',15); %black = 0% of patients had positive activity; red = 100% of patients had positive activity
        plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'o','color',cm_percent(percentPositive(v_m),:),'MarkerFaceColor',cm_percent(percentPositive(v_m),:),'MarkerSize',2.5); %black = 0% of patients had positive activity; red = 100% of patients had positive activity
        disp([num2str(round(v_m/nvert*100,2)) '%'])
    end
end
%SAVE FIGURE WITH COORDINATES
exportgraphics(gcf, [fig_name_vert_pos, '.png'])

%  BRAIN OF NUMBER OF PATIENTS NEAR EACH ELECTRODE
fig_name_vert_num = [hem, ' ', sx_input{1,1}, ' ', num2str(minnumpts), ' patients electrode coverage'];

if mp_count == 1
    figure('color','w','position',[230 171 1440 796],'Name',fig_name_vert_num);
    cm_npt = cbrewer2('Purples',7,'seq');
    
    
    getbrain4_ns('MNI','',1,0,hem,opscea_path,data_path) %plot empty brain
    shading interp
    axis equal; axis off; set(gca,'clipping','off'); lightsout; 
    litebrain(hem,.75)
    hold on
    
    for v_m = 1:nvert %for this specific vertex
    %    if nptAtVerts(v_m)>=minnumpts %if there are at least minnumpts number of patients with electrodes dst_radius away from the vertex
        if nptAtVerts(v_m)>=1 %if there are at least minnumpts number of patients with electrodes dst_radius away from the vertex
            plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'o','color', cm_npt(nptAtVerts(v_m),:),'MarkerFaceColor',cm_npt(nptAtVerts(v_m),:),'markersize',2.5); %pink lemonade plot of number of patients with vertices dst_radius away from electrodes
            disp([num2str(round(v_m/nvert*100,2)) '%'])
        end
    end
end

exportgraphics(gcf, [fig_name_vert_num, '.png'])


cd(opscea_path) %place so you don't have to change paths every time you run the code


