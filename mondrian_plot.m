function [sem_start,plot_start,plot_end] = mondrian_plot(pt,sz,perdur,ts_sx,to_plot,opscea_path,data_path)

% pt is a string such as 'UCSF4' or 'JaneDoe', acts as a prefix for files below
% sz is a string for '01' or other number of seizure for your patient, acts
% perdur is an integer indicating seconds before and after the analysis of
% the first symptom onset to analyze
% to_plot is whether or not to run the function

if to_plot
    % 1. Load data and find number of columns
    if ~exist(data_path,'dir'); error('Directory for your data needs to be corrected'); end

    cd([data_path pt '/' pt '_' sz]);

    sem_matrix_filename = [pt '_' sz '_mat.csv'];
    sem_matrix = readtable(sem_matrix_filename); %JK
    
    [~,t_cols] = size(sem_matrix); %check if transposed properly 
    
    % ----------------------------------------------------------------
    % 2. For loop to separate string features 
    y_label_names = {};      
    feature_el_vec = []; 
    sx_sem_matrix = [];
    sx_count = 0;


    first_vec = [];

    for i = 1:t_cols                                
        feature_el=sem_matrix.Properties.VariableNames{i}; %JK
        feature_el_vec = [feature_el_vec;sem_matrix.Properties.VariableNames{i}];
        anatomy = feature_el(1:2);
        position = feature_el(3);
%         motor = feature_el(4);
    
        % Anatomy
        switch anatomy 
            case 'lu'; full_anat = 'L Arm'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            case 'ru'; full_anat = 'R Arm'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            case 'lh'; full_anat = 'L Head'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            case 'rh'; full_anat = 'R Head'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            case 'le'; full_anat = 'L Eye '; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            case 're'; full_anat = 'R Eye'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            case 'lm'; full_anat = 'L Mouth'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            case 'rm'; full_anat = 'R Mouth'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
%             case 'ht'; full_anat = 'Head Turn'; sx_count = sx_count + 1; sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
%             case 'tt'; full_anat = 'Torso Turn';
%             case 'vx'; full_anat = 'Voice';
%             case 'gm'; full_anat = 'Gyratory Movement';
%             case 'rx'; full_anat = 'Rocking';
%             case 'bb'; full_anat = 'Bimanual Bipedal Automatism';
%             case 'wx'; full_anat = 'Walking';
%             case 'fx'; full_anat = 'Falling';
%             case 'px'; full_anat = 'Pedaling';
%             case 'ln'; full_anat = 'Left Hand';
%                 
%             case 'rn'; full_anat = 'Right Hand';
            case 'll'; full_anat = 'L Leg';
            case 'rl'; full_anat = 'R Leg';
%             case 'lf'; full_anat = 'Left Foot';
%             case 'rf'; full_anat = 'Right Foot';
%             case 'br'; full_anat = 'Behavioral Arrest';
%             case 'fa'; full_anat = 'Facial Automatism';
%             case 'oa'; full_anat = 'Oral Automatism';
%             case 'fa'; full_anat = 'Quadritonic';
%             case 'cg'; full_anat = 'Chapeau de Gendarme';
%             case 'fe'; full_anat = 'Facial Expression';
        end
    
        % Position
        switch position 
            case 'p'; full_pos = 'Proximal';
            case 'd'; full_pos = 'Distal';
%             case 'l'; full_pos = 'Left';
%             case 'r'; full_pos = 'Right';
%             case 'c'; full_pos = 'Center';
%             case 't'; full_pos = 'Twitch';
%             case 'y'; full_pos = 'Pull'; %y stands for yank so pull doesn't get confused with proximal
%             case 's'; full_pos = 'Superior';
%             case 'i'; full_pos = 'Inferior';
%             case 'f'; full_pos = 'Forward';
%             case 'b'; full_pos = 'Backward';
%             case 'n'; full_pos = 'Nonverbal';
%             case 'v'; full_pos = 'Verbal';
            case 'x'; full_pos = [];
        end
        y_label_names{sx_count} = [full_anat ' ' full_pos];

    end           

    % -----------------------------------------------------------------
    % 3. Separate numerical elements 
    nums_t_mat=sx_sem_matrix(:,:); % Separate numerical elements
    % -----------------------------------------------------------------
    % 5. Bin present symptoms (clean mat) and first appearance of
    % symptoms (ll_weight_times)

    [rows,cols] = size(nums_t_mat);

    t_sec = rows * .2; %convert to sec 
 
    col_num = [];

    first_auto = [];
    last_auto = [];

    first_tonic = [];
    last_tonic = [];

    first_clonic = [];
    last_clonic = [];

    cd(opscea_path)

    for n = 1:cols
        if ismember(1,nums_t_mat(:,n))
            clean_mat(:,n) = nums_t_mat(:,n); %collect all semiology with 1, 2, or 3 values in col vector
            %add first symptom start time
            first_auto(n) = find(clean_mat(:,n)==1,1,'first'); % Automatism
            last_auto(n) = find(clean_mat(:,n)==1,1,'last'); % Automatism
        end
        if ismember(2,nums_t_mat(:,n))
            clean_mat(:,n) = nums_t_mat(:,n); %collect all semiology with 1, 2, or 3 values in col vector

            first_tonic(n) =  find(clean_mat(:,n)==2,1,'first'); % Tonic
            last_tonic(n) = find(clean_mat(:,n)==2,1,'last'); % Tonic
        end
        if ismember(3,nums_t_mat(:,n))
            clean_mat(:,n) = nums_t_mat(:,n); %collect all semiology with 1, 2, or 3 values in col vector

            first_clonic(n) = find(clean_mat(:,n)==3,1,'first'); % Clonic
            last_clonic(n) = find(clean_mat(:,n)==3,1,'last'); % Clonic
        end
            %             first_vec = [first_vec min([first_auto,first_tonic,first_clonic])];
        if ~ismember(1,nums_t_mat(:,n)) && ~ismember(2,nums_t_mat(:,n)) && ~ismember(3,nums_t_mat(:,n))
            y_label_names{n}=''; %Delete y labels of row vectors without 1,2, or 3 value 
        end

    end

    bin_any = any(clean_mat);

    % -----------------------------------------------------------------
    % 6. Plot with ImageSC
    sem_start = round(min([first_auto(first_auto > 0) first_tonic(first_tonic>0) first_clonic(first_clonic>0)])*.2);
    sem_end = round(max([last_auto(last_auto > 0) last_tonic(last_tonic>0) last_clonic(last_clonic>0)])*.2);

    m=clean_mat(:,bin_any~=0)'; 
    semts=linspace(1,t_sec,rows);

    figure('Color','w','Name',['Semiology: ' pt '_' sz])
    imagesc(semts,1:size(m,1),m)

    hold on;
    plot_start = sem_start-(perdur*3);
    plot_end = sem_end+perdur;

    xlim([plot_start plot_end])

    colorbar('location','southoutside','Ticks',[.5:1:4.5],'TickLabels',{'No Motion','Automatism','Tonic','Clonic','Out of Video'}); %Place colorbar beneath graph; Hard code tick location by summing up # of values used (0-4 in this case)            
    cmap = [0.8,0.8,0.8;1,1,0;0,.4,1;.9,0,0;1,1,1]; % 0=Grey, 1=yellow, 2=blue, 3=red, 4=white
    colormap(gca,cmap);

    xlabel('Time (seconds)')

    y_ticks = 1:length(find(bin_any~=0));

    set(gca,'ytick',y_ticks,'CLimMode', 'manual', 'CLim', [0 5],'FontSize',12);
    set(gca,'YAxisLocation','right')


    y_labels = y_label_names(bin_any~=0);
    yticklabels('manual'); %remove 
    yticklabels(y_labels)
    

    set(gca,'YAxisLocation','right')

    xline([1.5:1:max(get(gca,'ytick'))],'k-',2.5)
    iter = 7.5; %for wide plot

    yline([min(xlim):iter:max(xlim)],'k-',.01)


    max(get(gca,'YTick'));

    yline(ts_sx,'m-',2) % vertical red line that follows time

end
