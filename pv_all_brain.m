function pv_all_brain(sx_input,lat_sxmx,npt,num_elecs,min_elec,minnumpts,opscea_path,data_path)

% created by Natalia Sucher 2/9/2023
% edited by Natalia Sucher 8/13/2023 

req_elecs = num_elecs >= min_elec; %matrix anatomies by patient with sufficient number of electrodes
[~,req_col] = size(req_elecs);

% establish neuroanatomy labels
fig_name_pv_all_pos = [sx_input{1,1} ': positive significance'];


figure('Name',fig_name_pv_all_pos,'Color','w'); % different figure for each symptom/mode combination
getbrain4_ns('MNI','',1,0,'r',opscea_path,data_path);
shading flat

npt_pos = nan(26,1);
req_pos = nan(27,npt);
percent_pos_req = nan(27,1);

cm_percent = cbrewer2('Reds',150,'cubic'); %color map 
cm_percent = cm_percent(20:120,:);
cm_npt = cbrewer2('Purples',7,'seq');

num_elecs_min = sum(req_elecs,2); % number of patients that have at least 5 electrodes in neuroanatomy
pos_num_elecs_min = num_elecs_min;
pos_num_elecs_min(pos_num_elecs_min < minnumpts) = NaN;

cd(data_path)

for l_i = 1:length(lat_sxmx)
    if l_i == 1
        pos_pv_T = readtable('pos_pv.xlsx', 'Sheet', lat_sxmx{l_i},'VariableNamingRule','preserve');
        pos_labels(:,1) = table2array(pos_pv_T(1:end,1));
        pos_pv_m = table2array(pos_pv_T(1:end,2:end));%(1:end,1:end)); % pvals per ll meandiff of neurosem across all electrodes
    end
    if req_col >= l_i
        req_idx = find(req_elecs(:,l_i) > 0);
        for r = 1:length(req_idx)
            if pos_pv_m(req_idx(r),l_i) < .05
                req_pos(req_idx(r),l_i) = pos_pv_m(req_idx(r),l_i);
            end
        end
    end
end



for label_i = 1:length(pos_labels)
    if sum(any(req_pos(label_i,:))) > 0 
        npt_pos(label_i) = sum(~isnan(req_pos(label_i,:)));
    end
end

npt_pos = npt_pos';
clear label_i

%PLOT RED BRAIN OF POSITIVE SIGNIFICANCE
cd(opscea_path)

for label_i = 1:length(npt_pos)
%     percent_pos(label_i) = npt_pos(label_i) / npt_pos(label_i);
    if pos_num_elecs_min(label_i) >= minnumpts
%             percent_pos_req(label_i) = npt_pos(label_i) / npt;
            percent_pos_req(label_i) = npt_pos(label_i) / pos_num_elecs_min(label_i);
            color_idx = round((npt_pos(label_i) / pos_num_elecs_min(label_i)) * length(cm_percent));
            if color_idx > 0
                highlightbrain('MNI',pos_labels(label_i),[1 1 1],[0 1],0,0,'r',opscea_path,data_path);
                highlightbrain('MNI',pos_labels(label_i),[cm_percent(color_idx,:);0 0 0],[0 1],0,0,'r',opscea_path,data_path);
            elseif color_idx == 0
                color_idx = 1;
                highlightbrain('MNI',pos_labels(label_i),[cm_percent(color_idx,:);0 0 0],[0 1],0,0,'r',opscea_path,data_path);
            end
    end
end

exportgraphics(gcf, [fig_name_pv_all_pos, '.png'])


%PLOT PURPLE BRAIN WITH # OF PATIENTS WITH MORE THAN 5 ELECTRODES
fig_name_pv_all_num_pts = [sx_input{1,1} ':number of patients with electrode coverage'];

figure('Name',fig_name_pv_all_num_pts,'Color','w'); % different figure for each symptom/mode combination
getbrain4_ns('MNI','',1,0,'r',opscea_path,data_path);
shading flat

for label_i = 1:length(npt_pos)
%     if num_elecs_min(label_i) >= minnumpts
    color_idx = round((num_elecs_min(label_i) / npt) * length(cm_npt));
    if color_idx > 0
        highlightbrain('MNI',pos_labels(label_i),[cm_npt(color_idx,:);0 0 0],[0 1],0,0,'r',opscea_path,data_path);
    end
%     end 
end

savefig([cd '/fig files/' fig_name_pv_all_num_pts])
exportgraphics(gcf, [cd '/png files/' fig_name_pv_all_num_pts,'.png'])
