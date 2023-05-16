function max_avg_MNI(sz_nns_mat,sz_w8s_mat,mni_xyz_cell,npt,hem,dst_radius)


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


mni_mesh_path = '/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/OPSCEADATA/avg_change_folders/cvs_avg35_inMNI152/Imaging/Meshes/';

cd('/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/OPSCEADATA/avg_change_folders/cvs_avg35_inMNI152/Imaging/Meshes/');


%START
load([mni_mesh_path 'cvs_avg35_inMNI152_' hem 'h_pial.mat']); %load right mni vertices
mni_h = cortex.vert;
nvert=size(mni_h,1);

cd('/Users/nataliasucher/Desktop/UCSF/Coding/OPSCEA')

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



minnumpts=1;
nptAtVerts=sum(haselecclose,2);
nptPositive=sum(max_w8s_rad>0,2); 
percentPositive=nptPositive./nptAtVerts; 


% plot3(mni_h(nptAtVerts<minnumpts,1),mni_h(nptAtVerts<minnumpts,2),mni_h(nptAtVerts<minnumpts,3),'.','color',[0 0 0],'markersize',3);
% figure('color','w','position',[230 171 1440 796]);
% getbrain4_ns('MNI','',1,0,hem) %plot empty brain
% shading interp
% axis equal; axis off; set(gca,'clipping','off'); lightsout; litebrain(hem,1)
% hold on
% 
% pause(1)

%RASPBERRY BRAIN OF % OF PATIENTS WITH POSITIVE ACTIVITY
% 
% for v_m = 1:nvert %for this specific vertex
%     if nptAtVerts(v_m)>=minnumpts %if there are at least minnumpts number of patients with electrodes dst_radius away from the vertex
% %         figure(f1)
% %         plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'.','color', cm(nptAtVerts(v_m),:),'markersize',15); %pink lemonade plot of number of patients with vertices dst_radius away from electrodes
%         plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'.','color',[percentPositive(v_m) 0 0],'markersize',15); %black = 0% of patients had positive activity; red = 100% of patients had positive activity
%         disp([num2str(round(v_m/nvert*100,2)) '%'])
%     end
% end
% 



%PINK LEMONADE BRAIN OF NUMBER OF PATIENTS NEAR EACH ELECTRODEE

cm = colormap(flipud(colormap(spring(npt))));


figure('color','w','position',[230 171 1440 796]);
getbrain4_ns('MNI','',1,0,hem) %plot empty brain
shading interp
axis equal; axis off; set(gca,'clipping','off'); lightsout; litebrain(hem,1)
hold on

for v_m = 1:nvert %for this specific vertex
    if nptAtVerts(v_m)>=minnumpts %if there are at least minnumpts number of patients with electrodes dst_radius away from the vertex
%         figure(f1)
        plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'.','color', cm(nptAtVerts(v_m),:),'markersize',15); %pink lemonade plot of number of patients with vertices dst_radius away from electrodes
%         plot3(mni_h(v_m,1),mni_h(v_m,2),mni_h(v_m,3),'.','color',[percentPositive(v_m) 0 0],'markersize',15); %black = 0% of patients had positive activity; red = 100% of patients had positive activity
        disp([num2str(round(v_m/nvert*100,2)) '%'])
    end
end

return


