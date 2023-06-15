% Natalia Sucher in the Kleen Lab, UCSF
% Created 1/31/2023
% Edited 5/25/2023


%
% Bug: cuts off neuroanatomy for lingual gyrus in activity_change.py or
% activity_plot.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE INPUTS DEPENDING ON YOUR OWN PREFERENCE


mni = 1; %average brain
% mni = 0; %individual brains

tic

% EDIT THIS TO REFLECT YOUR PATH

%   1. specify directory
opscea_path = '/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/';
data_path= [opscea_path 'OPSCEADATA/'];   %path for parameters sheet

 
% EDIT THIS TO REFLECT THE SYMPTOM
% sx_input = {'lhx','rhx','lud','rud', 'lup','rup'};
sx_input = {'chx'}; %{'cup'};; %{'chx'}; 
% sx_input = {'cex','cnx','cmx','cup','cud','cld','clp','fax','oax'};

% EDIT MODE (1 = AUTOMATISM, 2 = TONIC, 3 = CLONIC)
mx_input = {'2'};

% EDIT # OF SECONDS BEFORE AND AFTER SYMPTOM TO ANALYZE
perdur_input = '10';


% UNCOMMENT AND EDIT THIS IF YOU WANT TO SPECIFY PATIENT SEIZURES MANUALLY
manual_ptsz = {'EC91_03','EC96_01','EC107_01','EC133_03', 'EC166_01','EC228_03','EC229_02'}; % specify which patient and seizure

% MINIMUM ELECTRODES REQUIRED PER NEUROANATOMICAL LOCATIONS
min_elec = 5;

% MINIMUM & MAXIMUM PATIENTS REQUIRED FOR MAJORITY
min_pt = 4;
max_pt = 7;

% KEEP TRACK OF FOR LOOPS 
sz_count = 0;
sxmx_count = 0;

% CHOOSE NEUROANATOMICAL LOCATION TO ANALYZE (FREE SURFER LABELS)
reg_all = {'stg','mtg','itg','fus','tp','ph','hp','am','ent','pt','po','por','cmf','rmf','sf','sm','mof','lof','fp','prec','postc'};
% reg_all = {'prec'};

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(opscea_path)

% py_version('/Users/nataliasucher/opt/anaconda3/bin/python3.9');

%   1.5. for loop 
lat_sxmx = {};


elec_book = nan(200,3,length(manual_ptsz));
vert_book = nan(300000,3,length(manual_ptsz));

sz_nns_mat =  cell(1,length(manual_ptsz));
sz_w8s_mat = cell(1,length(manual_ptsz)); 
szxyz_mat = cell(1,length(manual_ptsz));

mni_xyz_cell = cell(1,length(manual_ptsz));

elecs_msize = zeros(2,length(manual_ptsz));
e_max_vec = ones(length(reg_all),length(manual_ptsz));

num_elecs = [];

for sx_i = 1:length(sx_input) % for loop throughout symptoms
    sx_name = sx_input{sx_i};

    for mx_i = 1:length(mx_input) % for loop throughout modes
        mx_name = mx_input{mx_i};

        sxmx_name = [sx_name ' ' mx_name];

        sxmx_count = sxmx_count + 1;

        sz_count = 0;
        for ptsz_i = 1:length(manual_ptsz) % for loop throughout seizures

            pt_sxmx_name=sxmx_name;
            
            split_manual_ptsz = split(manual_ptsz{ptsz_i},'_');
            pt_name = split_manual_ptsz{1};
            sz_name = split_manual_ptsz{2};
            ptsz_name = [pt_name '-' sz_name];

            cd(opscea_path)

            mondrian_plot(pt_name,sz_name,perdur_input)

            if pt_sxmx_name(1)=='c'
                sz_count = sz_count + 1;
                
                %%%%%%%%%%%%%
                % split_ptsz(manual_ptsz,ptsz_i)
                %%%%%%%%%%%%%
               
        

                cd(opscea_path)

                load([data_path pt_name '/Imaging/elecs/clinical_elecs_all.mat'],'elecmatrix','anatomy');

                em1 = elecmatrix;
                clear elecmatrix
                %%%%%%%%%%%%%
                [e_row,~] = size(em1);

                isR=nansum(em1(:,1))>0; isL=isR~=1;

                if isR
                    pt_sxmx_name(1) = 'l';
                    lat_sxmx{sz_count} = pt_sxmx_name; % add to cell of symptoms with laterality
%                     lat_brain_vec{sz_count} = 'r';
                elseif isL
                    pt_sxmx_name(1) = 'r';
                    lat_sxmx{sz_count} = pt_sxmx_name; % add to cell of symptoms with laterality

%                     lat_brain_vec{sz_count} = 'l';
                end
                %%%%%%%%%%%%


%                 if mni == 1
                if exist([data_path pt_name '/Imaging/elecs/clinical_TDT_elecs_all_warped.mat']) == 2
                    load([data_path pt_name '/Imaging/elecs/clinical_TDT_elecs_all_warped.mat'],'elecmatrix');
                    mni_xyz = elecmatrix;
                elseif exist([data_path pt_name '/Imaging/elecs/clinical_elecs_all_warped.mat']) == 2
                    load([data_path pt_name '/Imaging/elecs/clinical_elecs_all_warped.mat'],'elecmatrix');
                    mni_xyz = elecmatrix;
                else
                    [elecmatrix, ~] = make_clinical_elec_all_warped(pt_name);
                    mni_xyz = elecmatrix;
                    %todo: delete nans from elecmatrix to be same size as szxyz -- follow activity_change.py (retranslate it into matlab)
                end
                save('mni_elecmatrix.mat','mni_xyz')
                clear elecmatrix;
% 
%                 else
%                     mni_xyz = em1;
%                 end


                cd(opscea_path)
                [laterality, w8s_array, anat_array, good_mni] = pyrunfile("activity_change.py", ["laterality", "w8s_array", "anat_array","good_mni"], sxmx_input=pt_sxmx_name, ptsz_input=ptsz_name, perdur_input=perdur_input, opscea_path=opscea_path, data_path=data_path, sz_count=sz_count, sxmx_count=sxmx_count, ptsz_i=ptsz_i, min_elec=min_elec,e_row=e_row,mni_xyz=mni_xyz);
            
                good_mni_list = good_mni.tolist();
                good_mni_mat = nan(length(good_mni_list)-1,3);
                len_good_mni = length(good_mni_list)-1;
    

                [~, w8_cell, anat_cell, sz_w8s, sz_nns, szxyz, loaf] = activity_plot(string(laterality), w8s_array, anat_array, pt_sxmx_name, ptsz_name, data_path, opscea_path, ptsz_i, pt_name, sz_name, lat_sxmx, len_good_mni);      
           


            % TO DO: DELETE THIS IF/ELSE CONDITIONAL, SMUSH INTO ONE
            else
                sz_count = sz_count + 1;

                lat_sxmx{sz_count} = pt_sxmx_name;

                split_manual_ptsz = split(manual_ptsz{ptsz_i},'_');
                pt_name = split_manual_ptsz{1};
                sz_name = split_manual_ptsz{2};
                ptsz_name = [pt_name '-' sz_name];
                
                cd(opscea_path);

                [laterality, w8s_array, anat_array, good_mni] = pyrunfile("activity_change.py", ["laterality", "w8s_array", "anat_array","good_mni"], sxmx_input=pt_sxmx_name, ptsz_input=ptsz_name, perdur_input=perdur_input, opscea_path=opscea_path, data_path=data_path, sz_count=sz_count, sxmx_count=sxmx_count, ptsz_i=ptsz_i, min_elec=min_elec,e_row=e_row,mni_xyz=mni_xyz);
                
                % heatmap of electrical activity change during symptom onset
                [~, w8_cell, anat_cell, sz_w8s, sz_nns, szxyz, loaf] = activity_plot(string(laterality), w8s_array, anat_array, pt_sxmx_name, ptsz_name, data_path, opscea_path, ptsz_i, pt_name, sz_name, lat_sxmx, len_good_mni); 



            end
            
            for w_anat = 2:length(w8s_array)
                num_elecs(w_anat-1,ptsz_i) = length(w8_cell{1,w_anat-1});
            end
            
            for m = 1:len_good_mni
                for xyz = 1:3
                    good_mni_mat(m,xyz) = good_mni_list{m}{xyz};
                end
            end
            %ELECTRODES BY REGION PLOT
            reg_count = 0;
            for reg_i = 1:length(reg_all)

                reg_count = reg_count+1;
                cd(opscea_path)

                if mni == 1
                    [e_max] = pt_brain_elecs(manual_ptsz{ptsz_i},"MNI",reg_all{reg_i},lat_sxmx{sz_count}, mni_xyz, anatomy, reg_count);
                else
                    [e_max] = pt_brain_elecs(manual_ptsz{ptsz_i},pt_name,reg_all{reg_i},lat_sxmx{sz_count}, em1, anatomy, reg_count);
                end

                e_max_vec(reg_i,ptsz_i) = e_max;
            end

            mni_xyz_cell{ptsz_i} = good_mni_mat;
            sz_nns_mat{1,ptsz_i} =  sz_nns;
            sz_w8s_mat{1,ptsz_i} = sz_w8s; 
            szxyz_mat{1,ptsz_i} = szxyz;
            elecs_msize(1,ptsz_i) = max(e_max_vec(:,ptsz_i));
        end
      end 
end 

close all 
cd(opscea_path)

minnumpts=4;

bin_bilat %pixel plot of collapsed bilateral hemisphere 



minnumpts=4;

pv_all_brain(lat_sxmx,length(manual_ptsz),num_elecs,min_elec,minnumpts) %p value heatmap of combined total patients by neuroanatomical region

max_avg_MNI(sz_nns_mat,sz_w8s_mat,mni_xyz_cell,length(manual_ptsz),'r',dst_radius,minnumpts) %vertex heatmap on right hemisphere of brain

max_avg_MNI(sz_nns_mat,sz_w8s_mat,mni_xyz_cell,length(manual_ptsz),'l',dst_radius,minnumpts) %vertex heatmap on left hemisphere of brain







toc