function [w8_cell, sz_w8s, sz_nns, szxyz] = activity_plot(laterality, w8s_array, anat_array, pt_sxmx_name, ptsz_name, data_path, opscea_path, ptsz_i, pt_name, sz_name, len_good_mni)

cd(data_path);

%convert python xlsx generation to matrices, get rid of pos and neg 

all_T = readtable('all_sign_change.xlsx', 'Sheet', pt_sxmx_name,'VariableNamingRule','preserve');
pos_T = readtable('pos_sign_change.xlsx', 'Sheet', pt_sxmx_name,'VariableNamingRule','preserve');
neg_T = readtable('neg_sign_change.xlsx', 'Sheet', pt_sxmx_name,'VariableNamingRule','preserve');

all_pv_T = readtable('all_pv.xlsx', 'Sheet', pt_sxmx_name,'VariableNamingRule','preserve');
pos_pv_T = readtable('pos_pv.xlsx', 'Sheet', pt_sxmx_name,'VariableNamingRule','preserve');
neg_pv_T = readtable('neg_pv.xlsx', 'Sheet', pt_sxmx_name,'VariableNamingRule','preserve');

em_T = readtable('elec_matrix.xlsx','VariableNamingRule','preserve');

% all_m = table2array(all_T(~isnan(all_T(:,ptsz_i+1)),ptsz_i+1)); %extract linelength meandiff (positive or negative value of deviation from mean line length) 
all_m = table2array(all_T(1:end,ptsz_i+1)); %extract linelength meandiff (positive or negative value of deviation from mean line length) 
pos_m = table2array(pos_T(1:end,ptsz_i+1)); %extract linelength meandiff (positive value of deviation from mean line length) 
neg_m = table2array(neg_T(1:end,ptsz_i+1)); %extract linelength meandiff (negative value of deviation from mean line length) 

% all_m = table2array(all_T(~isnan(all_T(:,ptsz_i+1)),ptsz_i+1)); %extract linelength meandiff (positive or negative value of deviation from mean line length) 

pv_m = table2array(all_pv_T(1:end,ptsz_i+1)); % pvals per ll meandiff of neurosem across all electrodes
pos_pv_m = table2array(pos_pv_T(1:end,ptsz_i+1)); % pvals per ll meandiff of neurosem across all electrodes
neg_pv_m = table2array(neg_pv_T(1:end,ptsz_i+1)); % pvals per ll meandiff of neurosem across all electrodes

em_m = table2array(em_T(1:end,2:end));
sign_m = {all_m};


cd(opscea_path)

if strcmpi(laterality,pt_sxmx_name(1)) ~= 1 %contralateral only 

    class_wa = string(class(w8s_array));
    if class_wa ~= 'py.NoneType'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pv_brain(pos_pv_m,pos_pv_T,...
                 em_m,pt_name,sz_name,ptsz_name,ptsz_i,pt_sxmx_name,laterality,opscea_path,data_path);

        py_w8s_cell = cell(w8s_array);
        py_anat_cell = cell(anat_array);
        
        anat_cell = {};
        w8_cell = {};
    
        for w = 1:length(py_w8s_cell)
            for c = 1:length(py_w8s_cell{w})
                w8_cell{w}{c} = double(py_w8s_cell{w}{c});
                anat_cell{w} = string(py_anat_cell{w}{c});
            end
        end
    
    
        sign_now = sign_m{1,1};
        
        figure('Name',[pt_sxmx_name ' ' ptsz_name ': activity change'],'color','w'); % different figure for each symptom/mode combination
        subplot(2,1,1)

        cd(opscea_path);
        
        getbrain4_ns(pt_name,sz_name,1,0,laterality,opscea_path,data_path); %display brain 
        shading flat
        
        
        hold on;

        sz_w8s = squeeze(sign_now);
        sz_nns = ~isnan(sz_w8s);
        
        if any(pos_m)
            pos_mean_sz_LL = mean(pos_m(~isnan(pos_m(:,1))))*3.5;
        else
            pos_mean_sz_LL = mean(neg_m(~isnan(neg_m(:,1))))*3.5;
        end
        if any(neg_m)
            neg_mean_sz_LL = mean(neg_m(~isnan(neg_m(:,1))))*3.5;
        else
            neg_mean_sz_LL = mean(pos_m(~isnan(pos_m(:,1))))*3.5;
        end

        

        if exist("pos_mean_sz_LL","var") && exist("neg_mean_sz_LL","var")
            mean_sz_LL = mean([abs(pos_mean_sz_LL) abs(neg_mean_sz_LL)]);
        elseif exist("pos_mean_sz_LL","var") && ~exist("neg_mean_sz_LL","var")
            mean_sz_LL = pos_mean_sz_LL;
        elseif ~exist("pos_m","var") && exist("neg_mean_sz_LL","var")
            mean_sz_LL = abs(neg_mean_sz_LL);
        else
            mean_sz_LL = 20;
        end
        cax = [-mean_sz_LL mean_sz_LL];

        if ~isnan(mean_sz_LL)
        
        %FIND ELECMATRIX FOR SZ
            ptpath = [data_path pt_name];
            load([ptpath '/Imaging/elecs/clinical_elecs_all.mat'])
       
            %IMPORT PARAMS
            meshpath='/Imaging/Meshes/';
            Rcortex=load([ptpath meshpath pt_name '_rh_pial.mat']); 
                loaf.rpial=Rcortex; 
                Rcrtx=Rcortex.cortex; 
                clear Rcortex
            Lcortex=load([ptpath meshpath pt_name '_lh_pial.mat']); 
                loaf.lpial=Lcortex; 
                Lcrtx=Lcortex.cortex; 
                clear Lcortex
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~,prm_allPtSz]=xlsread([data_path 'OPSCEAparams.xlsx'],'params'); 
                fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
                prm=prm_allPtSz(strcmp(pt_name,prm_allPtSz(:,1))&strcmp(sz_name,prm_allPtSz(:,2)),:);
                if isempty(prm); error(['ATTENTION: No entry exists for ' pt_name ' seizure ' sz_name ' in the params master sheet']); end
        
            % Import parameters for patient's specific plot (layout of video frame)
            [~,plt]=xlsread([data_path 'OPSCEAparams.xlsx'],pt_name); 
                fields_PLOT=plt(1,:); plt(1,:)=[]; % header for columns of plotting parameters
                plottype=plt(:,strcmpi(fields_PLOT,'plottype')); %type of plot for each subplot (accepts: iceeg, surface, depth, or colorbar)
        
            cd 
        
            surfaces=plt(:,strcmpi(fields_PLOT,'surfaces'));
            surfacesopacity=plt(:,strcmpi(fields_PLOT,'surfacesopacity'));
        
            szxyz=[];

            idx=(ptsz_i-1)*3;
            
            szxyz(:,:)=em_m(1:len_good_mni,idx+1:idx+3);


            %% plot the weights on the brain
            params_gsp=str2double(prm{strcmp('gsp',fields_SZ)}); %gaussian spreading parameter (default 10)
            for j=1:size(plt,1) 
                      srf=regexp(surfaces{j},',','split'); % list the specific surfaces wanted for this subplot
                      srfalpha=regexp(surfacesopacity{j},',','split'); % list their corresponding opacities (values from 0 to 1; 0=invisible, 1=opaque)
                      if length(srf)~=length(srfalpha); msgbox('Number of surface to plot does not match number of alpha designations, check excel sheet'); return; end
                      acceptedterms={'rcortex','lcortex','rhipp','lhipp','ramyg','lamyg','wholebrain'};
                        for s=1:length(srf)
                            srf{s}=lower(srf{s}); %convert to lower case for easier string matching
                          if ~isempty(intersect(srf{s},acceptedterms)) %make sure user specified accepted terminologies for the meshes
                            switch srf{s} %see below for case "wholebrain"s
                                case 'rcortex'
                                    srfplot=Rcrtx; 
                                case 'lcortex'
                                    srfplot=Lcrtx; 
                            end
                          else
                            if laterality == 'r'
                                srfplot=Rcrtx;
                            elseif laterality == 'l'
                                srfplot=Lcrtx;
                            end
                          end 
                        end
                cd(opscea_path)
                if length(sz_nns) <= length(szxyz)
                    hh=ctmr_gauss_plot_edited(srfplot,szxyz(sz_nns,:),sz_w8s(sz_nns),cax,0,cmocean('balance'),params_gsp); 
                else 
%                     sz_w8s(1:length(szxyz)) = sz_w8s(sz_nns);
                    sz_w8s = sz_w8s(1:length(szxyz)); 
                    sz_nns = sz_nns(1:length(szxyz));
                    hh=ctmr_gauss_plot_edited(srfplot,szxyz(sz_nns(1:length(szxyz)),:),sz_w8s(sz_nns),cax,0,cmocean('balance'),params_gsp); 
                end
                colorbar("southoutside",'fontsize',18)
            end
            
            lightsout
            litebrain(laterality,1)
            hold on; 

    
            for i_row=1:size(szxyz(:,1))
                plot3(szxyz(i_row,1),szxyz(i_row,2),szxyz(i_row,3),'o','Color','k','MarkerFaceColor','k','MarkerSize',2.5); 
            end 
            

    %%%%%%%%%%%%%%
    % % anatomy and pval plots
        %to maintain y-labels with consistency across multiple patients
            

            reg_logic = zeros(1,length(anat_cell));
            mrkr = [];

            for u=1:length(pv_m)   
                if u < length(anat_cell)
                    if pos_pv_m(u)<.05
                        mrkr = [mrkr; 'ro']; 
                        reg_logic(u) = 1;
                    elseif isnan(pv_m(u))
                        reg_logic(u) = 0;
                    else 
                        mrkr = [mrkr; 'ko']; 
                        reg_logic(u) = 1;
                    end
                end
            end

            reg_idx = find(reg_logic); %index of neuroanatomy that has a significant p-value
            reg_plot = anat_cell(reg_idx); %string labels of neuroanatomical names
            w8s_plot = w8_cell(reg_idx);


            subplot(2,1,2)
            hold on;

            xlim(cax); %xlim([min(xlim)-.1 max(xlim)+.1]);


            for n=1:length(reg_plot)
                for w8 = 1:length(w8s_plot{n})       %plot individual electrodes (w8) per neurosemiology (u)
                    plot(w8s_plot{n}{w8},n*ones(length(w8s_plot{n}),1),mrkr(n,1:2)) % plot LL meandiff of each electrode
                end
                yline(n,':','LineWidth',.25,'Color',[.75,.75,.75]); % horizontal lines marking neuroanatomy
            end

            yticks(1:length(reg_plot))
            yticklabels(reg_plot)

            ylim([min(ylim)-1 max(ylim)+1]);

            xline(0,'k-'); % vertical line at x = 0 separating positive or negative activity

            set(gca,'FontSize',12)
            alpha(1)
        end
    else
        w8_cell = {};
        sz_nns = 0;
        szxyz=[];
        sz_w8s = [];

    end

end

%TRIM TRAILING NANS
if length(sz_nns) > length(szxyz)
    sz_nns = sz_nns(1:length(szxyz));
end

if length(sz_w8s) > length(szxyz)
    sz_w8s = sz_w8s(1:length(szxyz));
end

