function pv_all_brain(lat_sxmx,npt,num_elecs,min_elec,minnumpts)

%created by N. Sucher at the Kleen Lab in UCSF 2/9/2023


req_elecs = num_elecs > min_elec; %matrix anatomies by patient with sufficient number of electrodes
[~,req_col] = size(req_elecs);

% establish neuroanatomy labels

figure('Name',['positive significance'],'Color','w'); % different figure for each symptom/mode combination
getbrain4_ns('MNI','',1,0,'r');
shading flat

npt_pos = [];
req_pos = nan(27,length(lat_sxmx));
percent_pos = nan(27,1);
percent_pos_req = nan(27,1);

cm_percent = cbrewer2('Reds',150,'cubic'); %color map 
cm_percent = cm_percent(20:120,:);
% cm_percent = cm_percent(1:round((15*length(cm_percent))/16), :);
cm_npt = cbrewer2('Purples',7,'seq');

num_elecs_min = sum(req_elecs,2);


cd('/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/OPSCEADATA')

%assuming only one symptom is input
lat_count = 0;
for lat = 1:length(lat_sxmx)
    pos_pv_T = readtable('pos_pv.xlsx', 'Sheet', lat_sxmx{lat},'VariableNamingRule','preserve');
    pos_labels(:,1) = table2array(pos_pv_T(1:end,1));
%     [p_row,~] = size(pos_pv_T);
    pos_pv_m = table2array(pos_pv_T(1:end,2:end));%(1:end,1:end)); % pvals per ll meandiff of neurosem across all electrodes
    if req_col >= lat
        lat_count = lat_count + 1;
        req_idx = find(req_elecs(:,lat_count) > 0);
        for r = 1:length(req_idx)
            if pos_pv_m(req_idx(r),lat) < .05
                req_pos(req_idx(r),lat) = pos_pv_m(req_idx(r),lat);
            end
        end
    else
        k = 1;
    end
end

clear lat_count

% pos_pv_m = pos_pv_m(req_elecs);

for label_i = 1:length(pos_labels)
    if sum(any(req_pos(label_i,:))) > 0 
        npt_pos(label_i) = sum(~isnan(req_pos(label_i,:)));
    end
end

npt_pos = npt_pos';
clear label_i

%POSITIVE SIGNIFICANCE

cd('/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/')

for label_i = 1:length(npt_pos)
    percent_pos(label_i) = npt_pos(label_i) / npt_pos(label_i);
    if num_elecs_min(label_i) >= minnumpts
%             percent_pos_req(label_i) = npt_pos(label_i) / npt;
            percent_pos_req(label_i) = npt_pos(label_i) / num_elecs_min(label_i);
            color_idx = round((npt_pos(label_i) / num_elecs_min(label_i)) * length(cm_percent));
            if color_idx > 0
                highlightbrain('MNI',pos_labels(label_i),[1 1 1],[0 1],0,0,'r');
                highlightbrain('MNI',pos_labels(label_i),[cm_percent(color_idx,:);0 0 0],[0 1],0,0,'r');
            elseif color_idx == 0
                color_idx = 1;
                highlightbrain('MNI',pos_labels(label_i),[cm_percent(color_idx,:);0 0 0],[0 1],0,0,'r');
            end
    end
end
% colormap(cm_percent)
% cb1 = colorbar;
% cb1.TickLength = 0;

figure('Name',['number of patients'],'Color','w'); % different figure for each symptom/mode combination
getbrain4_ns('MNI','',1,0,'r');
shading flat

for label_i = 1:length(npt_pos)
%     if num_elecs_min(label_i) >= minnumpts
    color_idx = round((num_elecs_min(label_i) / npt) * length(cm_npt));
    if color_idx > 0
        highlightbrain('MNI',pos_labels(label_i),[cm_npt(color_idx,:);0 0 0],[0 1],0,0,'r');
    end
%     end 
end
% colormap(cm_npt)
% cb2 = colorbar;
% cb2.TickLength = 0;
% 

cd('/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/')

k = 1;