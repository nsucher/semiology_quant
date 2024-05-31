function [sem_start,plot_start,plot_end] = mondrian_plot(mx_input,sx_input,uber_ptsz,uber_lat,perdur,yes_plot,opscea_path,data_path,newfig)

% pt is a string such as 'UCSF4' or 'JaneDoe', acts as a prefix for files below
% sz is a string for '01' or other number of seizure for your patient, acts
% perdur is an integer indicating seconds before and after the analysis of
% the first symptom onset to analyze
% to_plot is whether or not to run the function

if yes_plot
    % 1. Load data and find number of columns
    if ~exist(data_path,'dir'); error('Directory for your data needs to be corrected'); end

    ptsz = split(uber_ptsz,'_');
    pt = ptsz{1};
    sz = ptsz{2};
    cd([data_path pt '/' pt '_' sz]);

    sem_matrix_filename = [pt '_' sz '_mat.csv'];
    sem_matrix = readtable(sem_matrix_filename); %JK
    
    [~,t_cols] = size(sem_matrix); %check if transposed properly 
    
    % ----------------------------------------------------------------
    % 2. For loop to separate string features 
    y_label_names = {};      
    % feature_el_vec = []; 
    sx_sem_matrix = [];
    sx_count = 0;


    for i = 1:t_cols
        cd('/scratch/semiology_quant');
        feature_el= sem_matrix.Properties.VariableNames{i}; %JK
        % feature_el_vec = [feature_el_vec;sem_matrix.Properties.VariableNames{i}];
        laterality = feature_el(1);
        anatomy = feature_el(1:2);
        position = feature_el(3);
        
        % Laterality
        if strcmpi(laterality,uber_lat{1}(1)) % plot only contralateral symptoms
        % Anatomy
            switch anatomy 
                case 'lu' 
                    full_anat = 'L Arm'; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'ru' 
                    full_anat = 'R Arm'; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'lh' 
                    full_anat = 'L Head'; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'rh' 
                    full_anat = 'R Head'; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'le' 
                    full_anat = 'L Eye '; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 're' 
                    full_anat = 'R Eye'; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'lm' 
                    full_anat = 'L Mouth'; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'rm' 
                    full_anat = 'R Mouth'; 
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'll' 
                    full_anat = 'L Leg';
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
                case 'rl' 
                    full_anat = 'R Leg';
                    sx_count = sx_count + 1; 
                    sx_sem_matrix(:,sx_count) = sem_matrix{:,i};
            end
        
    
            % Position
            switch position 
                case 'p' 
                    full_pos = 'Proximal';
                case 'd' 
                    full_pos = 'Distal';
                case 'x' 
                    full_pos = [];
            end

            y_label_names{sx_count} = [full_anat ' ' full_pos];

        end
        if strcmpi(sx_input{1}(2:3),feature_el(2:3)) % if the symptom is the one being analyzed
            sx_col = i;
            uber_lat1 = feature_el(1); % laterality of symptom onset that is being analyzed
        end

    end           

    % -----------------------------------------------------------------
    % 3. Separate numerical elements 
    % nums_t_mat_sx=sx_sem_matrix(:,:); % Separate numerical elements
    nums_t_mat=table2array(sem_matrix(:,:));
    % -----------------------------------------------------------------
    % 5. Bin present symptoms (clean mat) and first appearance of

    [rows,cols] = size(nums_t_mat);

    t_sec = rows * .2; %convert to sec 
 
    first_auto = [];
    sx_auto = [];
    last_auto = [];

    first_tonic = [];
    sx_tonic = [];
    last_tonic = [];

    first_clonic = [];
    sx_clonic = [];
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
        % if ~ismember(1,nums_t_mat(:,n)) && ~ismember(2,nums_t_mat(:,n)) && ~ismember(3,nums_t_mat(:,n))
        %     y_label_names{n}=''; %Delete y labels of row vectors without 1,2, or 3 value 
        % end

    end

    if ismember(1,nums_t_mat(:,sx_col))
        sx_auto = find(clean_mat(:,sx_col) == 1,1,'first');
    end

    if ismember(2,nums_t_mat(:,sx_col))
        sx_tonic = find(clean_mat(:,sx_col) == 2,1,'first');
    end

    if ismember(3,nums_t_mat(:,sx_col))
        sx_clonic = find(clean_mat(:,sx_col) == 3,1,'first');
    end


    bin_any = any(clean_mat);

    % -----------------------------------------------------------------
    % 6. Plot with ImageSC

    % start and stop from symptom onset and end
    
    sem_start = round(min([sx_auto(sx_auto > 0) sx_tonic(sx_tonic>0) sx_clonic(sx_clonic>0)])*.2);
    sem_end = round(max([last_auto(last_auto > 0) last_tonic(last_tonic>0) last_clonic(last_clonic>0)])*.2);


    plot_start = sem_start-(perdur*6);

    if plot_start < 1 
        plot_start = 1;
    elseif isempty(plot_start)
        plot_start = 1;
    end

    plot_end = sem_end;

    % x, y, and matrix values
    m=clean_mat(:,bin_any~=0)'; % matrix of seizure symptom mode over time
    semts=linspace(1,t_sec,rows); % time progression (x value)
    y = 1:size(m,1); % number of symptoms (y value)

    if newfig
        % Make figure and color
        sem_plot_name = [uber_lat{1} ' ' pt '_' sz ': semiology plot'];
        figure('Color','w','Name',sem_plot_name)
    end

    imagesc(semts,y,m)

    hold on;

    % x limits
    xlim([plot_start plot_end])

    % color
    cocolormaplorbar('location','southoutside','Ticks',.5:1:4.5,'TickLabels',{'No Motion','Automatism','Tonic','Clonic','Out of Video'}); %Place colorbar beneath graph; Hard code tick location by summing up # of values used (0-4 in this case)            
    cmap = [0.8,0.8,0.8;1,1,0;0,.4,1;.9,0,0;1,1,1]; % 0=Grey, 1=yellow, 2=blue, 3=red, 4=white
    colormap(gca,cmap);

    % y ticks
    y_ticks = 1:length(find(bin_any~=0));
    set(gca,'ytick',y_ticks,'CLimMode', 'manual', 'CLim', [0 5],'FontSize',12,'YAxisLocation','right');

    % y labesl
    y_labels = y_label_names(bin_any~=0);
    yticklabels(y_labels)
  
    % y line symptom separation
    for yl = 1.5:1:max(get(gca,'ytick'))
        line(xlim,[yl yl],'Color','k','LineWidth',2.5)
    end


    xlabel('Time (seconds)')

    % x line grids on and halfway between each x tick        
    
    x_ticks = (ceil(plot_start/10)*10):20:ceil(plot_end/10)*10;
    set(gca,'xtick',x_ticks)
    
    xt = xticks;
    iter = (xt(2)-xt(1))/2;
    if xt(length(xt)) < plot_end - iter
        xt(length(xt)+1) = xt(length(xt))+iter; % x line after final x tick
    end
    xt_iter = xt(1):iter:xt(length(xt)); % x line halfway between each x tick
    
    for xl = 1:length(xt_iter)
        line([xt_iter(xl) xt_iter(xl)],ylim,'Color','k','LineWidth',.001) % x line at each x tick
    end


end

if ~isempty(sem_start)
    line([sem_start sem_start],ylim,'Color','m','LineWidth',2.5)
end

% savefig([cd '/fig files/', sem_plot_name])
% exportgraphics(gcf, [cd '/png files/', sem_plot_name '.png'])
