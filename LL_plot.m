function [si,LL_s,ytl_LL,yt_LL,u2_s,u3_s]=LL_plot(new_anat,new_LL,ts,sem_start,plot_start,plot_end,sfx,data_path)
%     Created by Natalia Sucher and Jon Kleen May 10 2022, Updated May 26
%     2022 by NS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        %Electrode Activity
   
%writematrix(LL,'q8_LL.csv') %hardcoding specific file to feed to brain_w8s.m

%Electrode Activity   
  

           %Sort new anatomy
           [u1,u2,u3] = unique(new_anat(:,4)); % u1 = new_anat(u2); u2 = index of new_anat; u3 = index of u1
           [~,si] = sort(u3); % sort u3 by electrode index low to high
           u3_s = u3(si); % store sorted u3
           [~,u2_s] = unique(u3_s); %sort u2 (entries of first mentions) by deleting  repetitions in u3_s  
        
           %Sort linelength data 
           LL_s = new_LL(si,:); % sort LL data in ascending electrodes to plot

           abv_u1 = {};
           %Abbreviate ytick labels
           for k = 1:length(u1)
               %contains(u1{k},'superiortemporal')
               switch u1{k} 
                   case {'parstriangularis'}; abv_u1(k) = {'pt'};
                   case {'parsopercularis'}; abv_u1(k) = {'pop'};
                   case {'parsorbitalis'}; abv_u1(k) = {'por'};
                   %frontal
                   case {'rostralmiddlefrontal'}; abv_u1(k) = {'rmf'};
                   case {'caudalmiddlefrontal'}; abv_u1(k) = {'cmf'};
                   case {'lateralorbitofrontal'}; abv_u1(k) = {'lof'};
                   %case {'superiorfrontal'}; abv_u1(k) = {'sf'};
                   case {'medialorbitofrontal'}; abv_u1(k) = {'mof'};
                   %central
                   case {'precentral'}; abv_u1(k) = {'prec'};
                   case {'postcentral'}; abv_u1(k) = {'postc'}; 
                   %temporal
                   case {'middletemporal'}; abv_u1(k) = {'mtg'};
                   case {'superiortemporal'}; abv_u1(k) = {'stg'};                   
                   case {'inferiortemporal'}; abv_u1(k) = {'itg'};
                   case {'parahippocampal'}; abv_u1(k) = {'ph'};                   
                   case {'Right-Hippocampus'}; abv_u1(k) = {'rhp'};
                   case {'Left-Hippocampus'}; abv_u1(k) = {'lhp'};
                   case {'Right-Amygdala'}; abv_u1(k) = {'ram'};
                   case {'Left-Amygdala'}; abv_u1(k) = {'lam'};
                   case {'entorhinal'}; abv_u1(k) = {'ent'};
                   case {'fusiform'}; abv_u1(k) = {'fus'};                       
                   %poles    
                   case {'temporalpole'}; abv_u1(k) = {'tp'};
                   %case {'frontalpole'}; abv_u1(k) = {'fp'};
                   %marginal
                   case {'supramarginal'}; abv_u1(k) = {'sm'};
                   %other                       
                   case {'lingual'}; abv_u1(k) = {'lg'};
    %                case 'Right-Inf-Lat-Vent'; 
                   case {'Right-Cerebral-White-Matter'}; abv_u1(k) = {'rcwm'};
                   case {'Left-Cerebral-White-Matter'}; abv_u1(k) = {'lcwm'};
    %                case 'ctx...'; 
                   case {'bankssts'}; abv_u1(k) = {'bsts'};
                   case {'Right-choroid-plexus'}; abv_u1(k) = {'rchp'};
                   case {'Right-Putamen'}; abv_u1(k) = {'rput'};
                   case {'Right-VentralDC'}; abv_u1(k) = {'rvdc'};
                   case {'inferiorparietal'}; abv_u1(k) = {'ipt'};     
                   case {'superiorparietal'}; abv_u1(k) = {'spt'};                       

               end
           end
    
           %Average space between yticks for spacious labels
           [LL_s_row,~] = size(LL_s);
           [label_row,~] = size(u2_s);
           avg_yt_LL = nan(1,label_row);



           for i = 2:length(u2_s) + 1
               if i < length(u2_s)+1
                    avg_yt_LL(i-1) = u2_s(i-1) + (u2_s(i) - u2_s(i-1))/2; 
               elseif i == length(u2_s) + 1 %20
                   avg_yt_LL(i-1) = u2_s(i-1) + (LL_s_row - u2_s(i-1))/2;
               end
           end 



                               

           
           % variables for y axis 
           ytl_LL = abv_u1; %ytick labels for LL_s plot                    
           yt_LL = avg_yt_LL; %yticks for LL_s plot 
           
           % set new_LL and neuroanatomy labels as global
%            setGlobal_sem_w8s(LL_s,ytl_LL,yt_LL,u2_s,u3_s)

           % display plot with pcolor        
%            pcolor(ts(sem_start:end),1:size(new_LL,1),LL_s(:,sem_start:end));
           pcolor(ts(plot_start*sfx:plot_end*sfx),1:size(new_LL,1),LL_s(:,plot_start*sfx:plot_end*sfx));

%            cm = S.cm(floor(size(S.cm,1)/2):end,:);
           cd(data_path)
           cm = cmOPSCEAjet * .95;
           colormap(cm);
           shading flat;
          
           % modify axes
           xlabel('Time (seconds)')           
           set(gca,'ytick',yt_LL,'yticklabel',ytl_LL,'ydir','reverse','yaxislocation','right','fontsize',12)

           title('Line Length')
           tx=diff(xlim);
%            caxis(S.cax)
           caxis([0 20])

           hold on; plot(xlim,[u2_s u2_s],'w-')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        %Time Stamp
    % Time Marking
            %show vertical white lines to denote when semiology plot beginsand ends
%             SEMperiod = getGlobalSEMperiod;
%             plot([SEMperiod(1) SEMperiod(1)],ylim,'w-');
%             plot([SEMperiod(2) SEMperiod(2)],ylim,'w-');
            %hold on; 
            k = 2;

            