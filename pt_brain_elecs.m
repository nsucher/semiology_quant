function [e_max]=pt_brain_elecs(ptsz,pt_name,reg,pt_sxmx_name,elecmatrix,anatomy,reg_count,data_path)

%Edited by Natalia Sucher 8/13/23

% index of electrodes that are in the anatomical area (region) specified
% for the patient (pt) specified. 
% Note: em's 3 columns are ML(R+), AP(A+), DV (D+) 
% __common region names below:__
% 
% 'stg' --- 'superiortemporal'
% 'mtg' --- 'middletemporal'
% 'itg' --- 'inferiortemporal'
%
% 'ent' --- 'entorhinal'
% 'fus' --- 'fusiform'
% 'ph' --- 'parahippocampal'
% 'tp' --- 'temporalpole'
%
% 'hp' --- 'Right-Hippocampus' or
%          'Left-Hippocampus'
% 'am' --- 'Right-Amygdala' or
%          'Left-Amygdala'
%
% 'sm' --- 'supramarginal'
% 'prec' --- 'precentral'
% 'postc' --- 'postcentral'
% 'mof' --- 'medialorbitofrontal'
% 'lof' --- 'lateralorbitofrontal'
% 'pt' --- 'parstriangularis'
% 'pop' --- 'parsopercularis'
% 'por' --- 'parsorbitalis'
% 'rmf' --- 'rostralmiddlefrontal'
% 'sf' --- 'superiorfrontal'
% 'fp' --- 'frontalpole'

% pt_sxmx_name = symptom code + mode 


if nanmean(elecmatrix(:,1))>0; hem='Right'; else hem='Left'; end; hem=[hem '-']; 
isc=0;
if ischar(reg); reg={reg}; isc=1; end
for i=1:length(reg);
z=strcmp(reg{i},'stg'); if any(z); reg(i)={'superiortemporal'}; end
z=strcmp(reg{i},'mtg'); if any(z); reg(i)={'middletemporal'}; end
z=strcmp(reg{i},'itg'); if any(z); reg(i)={'inferiortemporal'}; end
z=strcmp(reg{i},'fus'); if any(z); reg(i)={'fusiform'}; end
z=strcmp(reg{i},'tp'); if any(z); reg(i)={'temporalpole'}; end
z=strcmp(reg{i},'ph'); if any(z); reg(i)={'parahippocampal'}; end
z=strcmp(reg{i},'hp'); if any(z); reg(i)={[hem 'Hippocampus']}; end
z=strcmp(reg{i},'am'); if any(z); reg(i)={[hem 'Amygdala']}; end
z=strcmp(reg{i},'ent'); if any(z); reg(i)={'entorhinal'}; end
z=strcmp(reg{i},'pt'); if any(z); reg(i)={'parstriangularis'}; end
z=strcmp(reg{i},{'po','pop'}); if any(z); reg(i)={'parsopercularis'}; end
z=strcmp(reg{i},{'por','porb'}); if any(z); reg(i)={'parsorbitalis'}; end
z=strcmp(reg{i},'cmf'); if any(z); reg(i)={'caudalmiddlefrontal'}; end
z=strcmp(reg{i},'rmf'); if any(z); reg(i)={'rostralmiddlefrontal'}; end
z=strcmp(reg{i},'sf'); if any(z); reg(i)={'superiorfrontal'}; end
z=strcmp(reg{i},'sm'); if any(z); reg(i)={'supramarginal'}; end
z=strcmp(reg{i},'mof'); if any(z); reg(i)={'medialorbitofrontal'}; end
z=strcmp(reg{i},'lof'); if any(z); reg(i)={'lateralorbitofrontal'}; end
z=strcmp(reg{i},'fp'); if any(z); reg(i)={'frontalpole'}; end
z=strcmp(reg{i},'prec'); if any(z); reg(i)={'precentral'}; end
z=strcmp(reg{i},'postc'); if any(z); reg(i)={'postcentral'}; end
end
z=strcmp(reg,'mfg'); if any(z); reg={'caudalmiddlefrontal','rostralmiddlefrontal'}; isc=0; end
z=strcmp(reg,'cing'); if any(z); reg={'caudalanteriorcingulate','isthmuscingulate','posteriorcingulate','rostralanteriorcingulate'}; isc=0; end


if isc; reg=reg{1}; end

if ischar(reg) 
%     [~,~,anatomy]=getelecs(pt,1); 
    elecs=find(strcmp(anatomy(:,4),reg))';
else
    for i=1:length(reg)
%         [~,~,anatomy]=getelecs(pt,1); 
        elecs=[elecs find(strcmp(anatomy(:,4),reg{i}))'];
    end
end


% 
% % [em,~,anatomy]=getelecs(pt,2);
% [elecs,em,~] = getregionelecs(pt,reg);


% Find laterality

isR=nansum(elecmatrix(:,1))>0; 
isL=isR~=1;

if isR
    lat = 'r';
elseif isL
    lat = 'l';
end

if reg_count == 1
    getbrain4_ns(pt_name,'',1,1,lat);
else
    getbrain4_ns(pt_name,'',0,0,lat);
end
shading flat;
alpha(.6)

hold on; 

cd(data_path);

all_T = readtable('all_sign_change.xlsx', 'Sheet', pt_sxmx_name,'VariableNamingRule','preserve');
all_a = table2array(all_T);

clear elecs_all_size 
 
elecs_all_size = ones(length(elecs),1);

%DO NOT DELETE-- DOT SIZE PLOT
%dot enlargement factor
dot_fact = 2;

for el = 1:length(elecs)
    if elecs(el) < length(all_a)
        t_el = table2array(all_T(elecs(el),{ptsz})).*dot_fact;
        if abs(t_el) < 1
            elecs_all_size(el,1) = 1;
        else
            elecs_all_size(el,1) = t_el;
        end
    end
end
elecs_all_size(isnan(elecs_all_size)) = 1;

for el = 1:length(elecs)
    if elecs_all_size(el,1) > 1
        mcolor = 'r'; % if electrode has positive activity change, make red
    elseif elecs_all_size(el,1) < 0
        mcolor = 'b'; % negative activity, make blue
    else
        mcolor = 'k'; % no activity, make black
    end
    h(el) = plot3(elecmatrix(elecs(el),1),elecmatrix(elecs(el),2),elecmatrix(elecs(el),3),'color', mcolor ,'marker','.','markersize', abs(elecs_all_size(el,1)));
end

e_max = max(abs(elecs_all_size));

if isempty(e_max)
    e_max = 1;
end



