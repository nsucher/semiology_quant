%uberfunction runs all code to generate figures for quantifying neural
%activity during seizure symptoms

%add mainpath variable (add to getbrain_ns and any ns or JON functions)

%INITIALIZE VARIABLES
% manual_ptsz = {'EC91_03','EC96_01','EC107_01','EC133_03', 'EC166_01','EC228_03','EC229_02'}; % specify which patient and seizure
manual_ptsz = {'EC229_02'}; % specify which patient and seizure

num_ptsz = length(manual_ptsz);

uber_pt = 'EC91'; %specific patient to show mondrian plot and opscea line length z-score graph for (in addition to all patients analyzed in neurosem_data)
uber_sz = '03'; %specific seizure

uber_ptsz = 'EC91_03';  
uber_i = find(strcmpi(manual_ptsz,uber_ptsz),1);

% index uber_pt/sz in manual_ptsz to match up with lat to only display contralateral side 

min_elec = 5; % MINIMUM ELECTRODES REQUIRED PER NEUROANATOMICAL LOCATIONS
min_pt = ceil(num_ptsz/2); %  MINIMUM PATIENTS REQUIRED FOR MAJORITY 

%SET PATH IF COPIED ELSEWHERE
uber_path = '/scratch/semiology_quant/';
cd(uber_path)
% cd '/data/shared/nsucher/share/semiology_quant/'


opscea_path = [pwd,'/'];
data_path= [opscea_path,'OPSCEADATA/'];   %path for parameters sheet
env_var = ["/home/",getenv("USER"),"/.conda/envs/sem_env/bin/python"];
env_path = join(env_var,""); %path for python environment

%DEBUG FOR NON NSUCHER 

% if contains(opscea_path,'/Users/nataliasucher/Desktop/UCSF/Coding/semiology_quant')
%     python_path = '/Users/nataliasucher/opt/anaconda3/bin/python';
% else
%     python_path = '/home/nsucher/.conda/envs/sem_quant3/bin/python';
% end


sx_input_list = {'chx'};% EDIT THIS TO REFLECT THE SYMPTOM
mx_input = {'2'};% EDIT MODE (1 = AUTOMATISM, 2 = TONIC, 3 = CLONIC)
perdur_input = 10; % EDIT # OF SECONDS BEFORE AND AFTER SYMPTOM TO ANALYZE


% APPROACH 1,2,3 AND ALL INDIVIDUAL PATIENT ACTIVITY CHANGES AND INDIVIDUAL AND INTRA-PATIENT PVALUES
for sx_i = 1:length(sx_input_list)
    sx_input = sx_input_list(sx_i);
    % [sx_plot,lat_sxmx] = neurosem_plot(uber_i,sx_input,mx_input,perdur_input,opscea_path,python_path,data_path,manual_ptsz,min_elec,min_pt,num_ptsz);
    [sx_plot,lat_sxmx] = neurosem_plot(sx_input,mx_input,perdur_input,opscea_path,data_path,manual_ptsz,min_elec,min_pt,num_ptsz,env_path);
    % yes_plot = 1; %1 = plot it, 0 = don't plot
    % for man_ptsz = 1:length(manual_ptsz)
    %     [sem_start,plot_start,plot_end] = mondrian_plot(sx_input,manual_ptsz{man_ptsz},lat_sxmx(man_ptsz),perdur_input,yes_plot,opscea_path,data_path);
    %     OPSCEA_sem_LL(manual_ptsz{man_ptsz},lat_sxmx(sx_i),1,sem_start,plot_start,plot_end,opscea_path,data_path) 
    % end
end
% uber_lat = lat_sxmx{uber_i}(1); %laterality of focused seizure

% SYMPTOM INDEX AS TIME SERIES (MONDRIAN PLOT)
% yes_plot = 1; %1 = plot it, 0 = don't plot
% 
% [sem_start,plot_start,plot_end] = mondrian_plot(uber_ptsz,uber_lat,perdur_input,yes_plot,opscea_path,data_path);
% 
% % OPSCEA ICEEG, LINE LENGTH TRANSFORM AS PLOT AND 3D BRAIN HEATMAP
% OPSCEA_sem_LL(uber_ptsz,uber_lat,1,sx_plot,plot_start,plot_end,opscea_path,data_path) 