% bin pseudocode

%   align and sectionalize left and right brain using grid of 1cm squares
%   bin majority percent in each square 
%   add up majority percent in left and right brains
%   create new  brain with summed majority percent in each grid as colors
%   

Y=-100:10:100;
Z=-100:10:100;
M=nan(length(Y),length(Z),7);


for p=1:num_ptsz
    em = szxyz_mat{p};
    w8s = sz_w8s_mat{p};
    if ~isempty(em)
        for y=1:length(Y)-1
            for z=1:length(Z)-1
                        insquare=find( em(:,2) > Y(y)   & ...
                                       em(:,2) < Y(y+1) & ...
                                       em(:,3) > Z(z)   & ...
                                       em(:,3) < Z(z+1));
                        if ~isempty(insquare)
                            if z > 1
                                M(y+1,z-1,p)=mean(w8s(insquare));
                            end
                        end
            end
        end
    else
        M(1,1,p) = nan;
    end
end





Mbin=M>5; %any ll mean difference (weight change) values above 5 for each patient averaged in each square
Mnpt=sum(~isnan(M),3);
Mnpt_req = sum(~isnan(M),3);
Mnpt_req(find(Mnpt < minnumpts)) = 0;
Mptsig=sum(Mbin,3);
Mptsig_req = sum(Mbin,3);
Mptsig_req(find(Mptsig < minnumpts)) = 0;
Mpercent = Mptsig_req'./Mnpt_req'*100;
Mpercent(Mpercent == 0) = 2;
Mpercent(isnan(Mpercent)) = 0;


figure('name',num2str(minnumpts),'color','w','position',[230 171 1000 796]);

% PLOT GRID FOR GLASS BRAIN
ylim(Y([1 end])); zlim(Z([1 end])); %in same axis limits for orientation
axis on
set(gca,'ytick',Y,'ztick',Z,'LineWidth',1,'GridColor','k') %grid lines where the boundaries should be

% PLOT 3D GLASS BRAIN
ax1 = subplot(2,2,1); 
getbrain4_ns('MNI','',1,0,'r',opscea_path,data_path); %brain for orientation
hold on;
shading interp
axis equal; axis off; set(gca,'clipping','off'); lightsout; litebrain('r',1)
alpha .5 % transparency to see the grid lines better

ax3 = subplot(2,2,3); % total significant patients per square
pcolor(Mptsig'); 
axis equal; axis off; 
% colormap(ax3,raspberry)
% figure('color','w')
% cm_sig = cbrewer2('PuRd',7,'seq');

cm_sig = cbrewer2_ns('PuRd',7,'seq');


cm_sig = [.85 .85 .85; cm_sig];
colormap(ax3,cm_sig)
% colorbar
pos_pt_data = flipud(Mptsig');


ax4 = subplot(2,2,4);  %percent sigificant patients per square
pcolor(Mpercent); 
axis equal; axis off;
cm_percent = cbrewer2_ns('Reds',150,'cubic'); %color map 
cm_percent = cm_percent(20:120, :);
cm_percent = [.85 .85 .85; cm_percent];

percent_pt_data = flipud(Mpercent);




% cm_percent = [.85 .85 .85; cm_percent];

%color map
colormap(ax4,cm_percent)

% colormap(ax4,pink_lemonade)

ax2 = subplot(2,2,2); %total patients per square
pcolor(Mnpt'); 
axis equal; axis off; 
% colormap(ax2,mintyfresh)
cm_npt = cbrewer2_ns('Purples',7,'seq');
cm_npt = [.85 .85 .85; cm_npt];
colormap(ax2,cm_npt)

num_pt_data = flipud(Mnpt');

k = 1;