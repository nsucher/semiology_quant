%uberfunction runs all code to generate figures for quantifying neural
%activity during seizure symptoms

uber_pt = 'EC91'; %specific patient to show mondrian plot and opscea line length z-score graph for (in addition to all patients analyzed in neurosem_data)
uber_sz = '03'; %specific seizure

[ts_sx, sx_sec] = neurosem_data(uber_pt);

[sem_start,plot_start,plot_end] = mondrian_plot(uber_pt,uber_sz,10,ts_sx/5,1);
% neurosem_data
 
OPSCEA_sem_LL(uber_pt,uber_sz,1,sem_start,ts_sx/5,plot_start,plot_end) 
