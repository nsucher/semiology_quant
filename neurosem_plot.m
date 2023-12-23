function [ts_sx,sx_sec] = neurosem_plot(uber_pt,sx_input,mx_input,perdur_input,opscea_path,python_path,data_path,manual_ptsz,min_elec,min_pt,num_ptsz)


% Function to generate figures analyzing change of neural activity during onset of seizure symptom
% output: ts_sx is the timestamp of the first onset of the seizure symptom


% Natalia Sucher in the Kleen Lab, UCSF
% Created 1/31/2023
% Edited 12/12/2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE INPUTS DEPENDING ON YOUR OWN PREFERENCE

mni = 1; %average brain
% mni = 0; %individual brains

tic
 

ts_sx_vec = []; %collect start time of symptom onset in seconds


% KEEP TRACK OF FOR LOOPS 
sxmx_count = 0;

% CHOOSE NEUROANATOMICAL LOCATION TO ANALYZE (FREE SURFER LABELS)
reg_all = {'stg','mtg','itg','fus','tp','ph','hp','am','ent','pt','po','por','cmf','rmf','sf','sm','mof','lof','fp','prec','postc'};

neuroanat_list = ['frontalpole', 'parstriangularis', 'parsopercularis', 'parstriangularis',...
    'parsopercularis', 'parsorbitalis', 'rostralmiddlefrontal', 'caudalmiddlefrontal',...
    'lateralorbitofrontal','superiorfrontal','medialorbitofrontal','precentral','postcentral',...
    'inferiorparietal','superiorparietal','supramarginal','temporalpole','middletemporal',...
    'superiortemporal','inferiortemporal','parahippocampal','Right-Hippocampus','Left-Hippocampus',...
    'Right-Amygdala','Left-Amygdala','entorhinal','bankssts','fusiform', 'lingual'];

abv_neuroanat_list = ['front-pole','parstri','parsop','parsorb','rost-midfront','caud-midfront',...
    'latorb-front','sup-front','medorb-front','precentral','postcentral','inf-par','sup-par',...
    'supra-marg','temp-pole','mid-temp','sup-temp','inf-temp','parahip','R-Hip','L-Hip','R-Amyg',...
    'L-Amyg','entorhinal','bankssts','fusiform','lingual'];

num_anat = length(neuroanat_list);

% CHOOSE RADIUS HERE
dst_radius = 10; % minimum distance in mm from electrode to each vertex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DELETE PREVIOUS XLSX FILES
cd(data_path)

delete pos_pv.xlsx
delete neg_pv.xlsx
delete all_pv.xlsx

delete all_sign_change.xlsx
delete pos_sign_change.xlsx
delete neg_sign_change.xlsx
delete cut_sign_change.xlsx

delete elec_matrix.xlsx
delete elec_w8s.xlsx

%%%%%%%%%%%%%%‰%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(opscea_path)

%   1.5. for loop 
lat_sxmx = {};

sz_nns_mat =  cell(1,num_ptsz);
sz_w8s_mat = cell(1,num_ptsz); 
szxyz_mat = cell(1,num_ptsz);

mni_xyz_cell = cell(1,num_ptsz);

elecs_msize = zeros(2,num_ptsz);
e_max_vec = ones(length(reg_all),num_ptsz);

num_elecs = [];
first_sx_vec = zeros(1,num_ptsz); %initialize empty vector to collect timestamps for onset of symtpoms for each patient

% for sx_i = 1:length(sx_input) % for loop throughout symptoms
%     sx_name = sx_input{sx_i};
sx_name = sx_input{1};
for mx_i = 1:length(mx_input) % for loop throughout modes
    mx_name = mx_input{mx_i};

    sxmx_name = [sx_name ' ' mx_name];

    sxmx_count = sxmx_count + 1;

    sz_count = 0;
    for ptsz_i = 1:num_ptsz % for loop throughout seizures

        pt_sxmx_name=sxmx_name;
        
        split_manual_ptsz = split(manual_ptsz{ptsz_i},'_');
        pt_name = split_manual_ptsz{1};
        sz_name = split_manual_ptsz{2};
        ptsz_name = [pt_name '-' sz_name];

        cd(opscea_path)

        sz_count = sz_count + 1;

        load([data_path pt_name '/Imaging/elecs/clinical_elecs_all.mat'],'elecmatrix','anatomy');

        em1 = elecmatrix;
        clear elecmatrix
        %%%%%%%%%%%%%
        [e_row,~] = size(em1);

        isR=nansum(em1(:,1))>0; isL=isR~=1;

        if isR
            pt_sxmx_name(1) = 'l';
            lat_sxmx{sz_count} = pt_sxmx_name; % add to cell of symptoms with laterality
        elseif isL
            pt_sxmx_name(1) = 'r';
            lat_sxmx{sz_count} = pt_sxmx_name; % add to cell of symptoms with laterality

        end
        %%%%%%%%%%%%


        if exist([data_path pt_name '/Imaging/elecs/clinical_TDT_elecs_all_warped.mat'],"file")
            load([data_path pt_name '/Imaging/elecs/clinical_TDT_elecs_all_warped.mat'],'elecmatrix');
            mni_xyz = elecmatrix;
        elseif exist([data_path pt_name '/Imaging/elecs/clinical_elecs_all_warped.mat'],"file") 
            load([data_path pt_name '/Imaging/elecs/clinical_elecs_all_warped.mat'],'elecmatrix');
            mni_xyz = elecmatrix;
        else
            [elecmatrix, ~] = make_clinical_elec_all_warped(pt_name,opscea_path,data_path);
            mni_xyz = elecmatrix;
        end
        save('mni_elecmatrix.mat','mni_xyz')
        clear elecmatrix;

        cd(opscea_path)

        terminate(pyenv)
        pyenv("ExecutionMode","OutOfProcess","Version",python_path);

        
        [laterality, w8s_array, anat_array, good_mni, first_mx] = pyrunfile("activity_change.py", ["laterality", "w8s_array", "anat_array","good_mni","first_mx"], sxmx_input=pt_sxmx_name, ptsz_input=ptsz_name, perdur_input=perdur_input, opscea_path=opscea_path, data_path=data_path, sz_count=sz_count, sxmx_count=sxmx_count, ptsz_i=ptsz_i, min_elec=min_elec,e_row=e_row,mni_xyz=mni_xyz);
        
%             all_pv_m = nan(num_anat,num_ptsz);%
%             all_sign_m = nan(num_anat,num_ptsz);%
%             elec_w8s_m = nan(1,num_ptsz); %will grow later--depends on patient with most electrodes
% 
% 
%             [laterality, w8s_array, anat_array, good_mni, first_mx] = plot_change(num_ptsz,pt_sxmx_name,sx_name,mx_name,ptsz_name,pt_name,perdur_input,opscea_path,data_path,sz_count,sxmx_count,ptsz_i,min_elec,e_row,mni_xyz,neuroanat_list,abv_neuroanat_list,num_anat,all_pv_m,all_sign_m,elec_w8s_m);

        if matches(pt_name,uber_pt)
            first_sx_vec(1,ptsz_i) = first_mx;
            ts_sx = int64(first_mx);
        end
      
%         ts_sx_vec(sx_i,ptsz_i) = ts_sx; %pt = cols; symptoms = rows
        good_mni_list = good_mni.tolist();
        good_mni_mat = nan(length(good_mni_list)-1,3);
        len_good_mni = length(good_mni_list)-1;


        [w8_cell, sz_w8s, sz_nns, szxyz] = activity_plot(string(laterality), w8s_array, anat_array, pt_sxmx_name, ptsz_name, data_path, opscea_path, ptsz_i, pt_name, sz_name, lat_sxmx, len_good_mni);      
       
        
        for w_anat = 2:length(w8s_array)
            num_elecs(w_anat-1,ptsz_i) = length(w8_cell{1,w_anat-1});
        end
        
        for m = 1:len_good_mni
            for xyz = 1:3
                good_mni_mat(m,xyz) = good_mni_list{m}{xyz};
            end
        end
        %ELECTRODES BY REGION PLOT
        % reg_count = 0;
        % for reg_i = 1:length(reg_all)
        % 
        %     reg_count = reg_count+1;
        %     cd(opscea_path)

%            if mni == 1
%                [e_max] = pt_brain_elecs(manual_ptsz{ptsz_i},"MNI",reg_all{reg_i},lat_sxmx{sz_count}, mni_xyz, anatomy, reg_count, opscea_path, data_path);
%            else
%                [e_max] = pt_brain_elecs(manual_ptsz{ptsz_i},pt_name,reg_all{reg_i},lat_sxmx{sz_count}, em1, anatomy, reg_count, opscea_path, data_path);
%            end

%            e_max_vec(reg_i,ptsz_i) = e_max;
%        end

        mni_xyz_cell{ptsz_i} = good_mni_mat;
        sz_nns_mat{1,ptsz_i} =  sz_nns;
        sz_w8s_mat{1,ptsz_i} = sz_w8s; 
        szxyz_mat{1,ptsz_i} = szxyz;
%        elecs_msize(1,ptsz_i) = max(e_max_vec(:,ptsz_i));
%         cd(opscea_path)
%         mondrian_plot(pt_name,sz_name,10,ts_sx/5,1,opscea_path,data_path);
    end
end 
 

close all 

mp_count = 0;

for minnumpts = 1:min_pt
    mp_count = mp_count + 1;
    cd(opscea_path)
    
    % PAPER: FIGURE __ / POSTER: APPROACH 3
    bin_bilat %pixel plot of collapsed bilateral hemisphere 

    % PAPER: FIGURE __/ POSTER: APPROACH 2 FIGURES VERTEX BY VERTEX
    % max_avg_MNI(sz_nns_mat,sz_w8s_mat,mni_xyz_cell,num_ptsz,'r',dst_radius,minnumpts,opscea_path,data_path,mp_count) %vertex heatmap on right hemisphere of brain
    % max_avg_MNI(sz_nns_mat,sz_w8s_mat,mni_xyz_cell,num_ptsz,'l',dst_radius,minnumpts,opscea_path,data_path,mp_count) %vertex heatmap on left hemisphere of brain

end
delete mp_count

cd(opscea_path)

pv_all_brain(sx_input,lat_sxmx,num_ptsz,num_elecs,min_elec,min_pt,opscea_path,data_path) %p value heatmap of combined total patients by neuroanatomical region

% [sem_start,plot_start,plot_end] = mondrian_plot(uber_pt,uber_sz,10,ts_sx/5,1,opscea_path,data_path);

%OPSCEA_sem_LL(uber_pt,uber_sz,1,sem_start,ts_sx/5,plot_start,plot_end) 
sx_sec = first_sx_vec/5;

% 
k=2;