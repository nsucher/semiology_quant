function pv_brain(pos_pv_m,neg_pv_m,pos_pv_T,neg_pv_T,em_m,pt_name,sz_name,ptsz_name,ptsz_i,sxmx_name,laterality)

%created by N. Sucher at the Kleen Lab in UCSF 2/9/2023
%edited 8/8/23

%         establishing brain pval color
whiterows=101;
max_pval=5;
lowestcolor=.85;

% establish neuroanatomy labels
pos_labels = table2array(pos_pv_T(1:end,1));
% neg_labels = table2array(neg_pv_T(1:end,1));

szxyz=[];

idx=(ptsz_i-1)*3;

szxyz(:,:)=em_m(:,idx+1:idx+3);

if any(pos_pv_m)
    figure('Name',[sxmx_name ' ' ptsz_name ': positive significance'],'Color','w'); % different figure for each symptom/mode combination


    getbrain4_ns(pt_name,sz_name,1,0,laterality);
   
    caxis([-max_pval max_pval])

    shading flat

    cm=makecm([lowestcolor lowestcolor lowestcolor;1 .75 .15],whiterows-1); %yellow for positive
    cm=[ones(whiterows,3); cm]; 

    for label_i = 1:length(pos_labels)
        pv_now = pos_pv_m(label_i);
        if pv_now < .05

            pv = -log10(pv_now);
            pv = pv/max_pval; %normalize pvalues to maxp

            pv(pv>1)=1; %any pvalues stronger than maxp reduced to maxp for plotting/colormap entry purposes
            pv=round(pv *(whiterows-1)); %get the row index for the main colormap
            pv=pv+whiterows; %adjust the row index for the filler whiterows in the colormap
            if ~isnan(pv)
                if pv == 101
                    pv = 102;
                end
                c = cm(pv,:);
                highlightbrain(pt_name,pos_labels(label_i),[c;0 0 0],[0 1],0,0,laterality);
                colormap(gca,cm)

                for i_row=1:size(szxyz(:,1))
                    plot3(szxyz(i_row,1),szxyz(i_row,2),szxyz(i_row,3),'k.','markersize',15); 
                end 
            end
        end
    end
end

% 
% if any(neg_pv_m)
% %     figure('Name',[sxmx_name ' ' ptsz_name ': negative significance'],'Color','w'); % different figure for each symptom/mode combination
% % 
% % 
% %     getbrain4_ns(pt_name,sz_name,1,0,laterality);
%    
% %     caxis([-max_pval max_pval])
% 
% %     shading flat
% % 
% %     cm=makecm([1 .75 .15;lowestcolor lowestcolor lowestcolor],whiterows-1); %yellow for positive
% %     cm=flipud(cm);
% %     cm=[ones(whiterows,3); cm]; 
% 
%     for label_i = 1:length(neg_labels)
%         pv_now = neg_pv_m(label_i);
%         if pv_now < .05
% 
%             pv = -log10(pv_now);
%             pv = pv/max_pval; %normalize pvalues to maxp
% 
%             pv(pv>1)=1; %any pvalues stronger than maxp reduced to maxp for plotting/colormap entry purposes
%             pv=round(pv *(whiterows-1)); %get the row index for the main colormap
%             pv=pv+whiterows; %adjust the row index for the filler whiterows in the colormap
%             if ~isnan(pv)
%                 if pv == 101
%                     pv = 102;
%                 end
%                 c = cm(pv,:);
% 
%                 highlightbrain(pt_name,neg_labels(label_i),[c;0 0 0],[0 1],0,0,laterality);
% 
%                 colormap(gca,cm)
% %                     log_pos_pv(label_i)=pv_now;
% 
%                 for i_row=1:size(szxyz(:,1))
%                     plot3(szxyz(i_row,1),szxyz(i_row,2),szxyz(i_row,3),'k.','markersize',15); 
%                 end 
%             end
%         end
%     end
% end
% 

