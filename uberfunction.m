%uberfunction runs all code to generate figures for quantifying neural
%activity during seizure symptoms

%add mainpath variable (add to getbrain_ns and any ns or JON functions)

%INITIALIZE VARIABLES
manual_ptsz = {'EC91_03','EC96_01','EC107_01','EC133_03', 'EC166_01','EC228_03','EC229_02'}; % specify which patient and seizure
num_ptsz = length(manual_ptsz);

uber_pt = 'EC91'; %specific patient to show mondrian plot and opscea line length z-score graph for (in addition to all patients analyzed in neurosem_data)
uber_sz = '03'; %specific seizure

min_elec = 5; % MINIMUM ELECTRODES REQUIRED PER NEUROANATOMICAL LOCATIONS

min_pt = ceil(num_ptsz/2); %  MINIMUM PATIENTS REQUIRED FOR MAJORITY 

%SET PATH IF COPIED ELSEWHERE
opscea_path = [pwd, '/'];
% opscea_path = '/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'; %set your own directory here
% opscea_path = '/Users/nataliasucher/Desktop/UCSF/Coding/semiology_quant';
python_path = '/home/nsucher/.conda/envs/semiology_quant/bin/python';
data_path= [opscea_path 'OPSCEADATA/'];   %path for parameters sheet
sx_input_list = {'chx'};% EDIT THIS TO REFLECT THE SYMPTOM
mx_input = {'2'};% EDIT MODE (1 = AUTOMATISM, 2 = TONIC, 3 = CLONIC)
perdur_input = '10'; % EDIT # OF SECONDS BEFORE AND AFTER SYMPTOM TO ANALYZE

for sx_i = 1:length(sx_input_list)
    sx_input = sx_input_list(sx_i);
    [ts_sx, sx_sec] = neurosem_plot(uber_pt,sx_input,mx_input,perdur_input,opscea_path,python_path,data_path,manual_ptsz,min_elec,min_pt,num_ptsz);
end
% 
% manual_ptsz = {'EC91_03','EC96_01','EC107_01','EC133_03', 'EC166_01','EC228_03','EC229_02'}; % specify which patient and seizure
% % 
% for m_p = length(manual_ptsz)
[sem_start,plot_start,plot_end] = mondrian_plot(uber_pt,uber_sz,10,ts_sx/5,1,opscea_path,data_path);
% end
 
OPSCEA_sem_LL(uber_pt,uber_sz,1,sem_start,ts_sx/5,110,250,opscea_path,data_path) 