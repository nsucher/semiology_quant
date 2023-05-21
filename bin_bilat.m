% bin pseudocode

%   align and sectionalize left and right brain using grid of 1cm squares
%   bin majority percent in each square 
%   add up majority percent in left and right brains
%   create new  brain with summed majority percent in each grid as colors
%   

%em must be the physical image? or the ll mean diff weights?



Y=-100:10:100;
Z=-100:10:100;
M=nan(length(Y),length(Z),7);


for p=1:7
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

Mbin=M>0; %any average weight values above zero for each patient in each square
Mnpt=sum(~isnan(M),3);
Mptsig=sum(Mbin,3);
Mpercent = Mptsig'./Mnpt'*100;



figure('color','w','position',[230 171 1000 796]);

ax1 = subplot(2,2,1); 
getbrain4_ns('MNI','',1,0,'r'); %brain for orientation
hold on;
shading interp
axis equal; axis off; set(gca,'clipping','off'); lightsout; litebrain('r',1)
% alpha(.6)

ylim(Y([1 end])); zlim(Z([1 end])); %in same axis limits for orientation
axis on
set(gca,'ytick',Y,'ztick',Z) %grid lines where the boundaries should be
alpha .5 % transparency to see the grid lines better

% colormap(ax1,colormap([1 1 1]))


% axis on
% set(gca,'ytick',Y,'ztick',Z) %grid lines where the boundaries should be
% grid on;
% axis on
% set(gca,'ytick',Y,'ztick',Z) %grid lines where the boundaries should be

%colors
% pink_lemonade = colormap(spring);
% pink_lemonade(1,:) = [1 1 1]; %set zero equal to white

% raspberry=makecm([0 0 0;0 1 1],7); %red-black colormap
% raspberry = [1 1 1; raspberry];

% 
% mintyfresh=makecm([0 1 1;0 0 0],7); %cyan-black colormap
% mintyfresh(1,:) = [1 1 1];
% 

% cm_npt = cbrewer2('Purples',7,'seq');
% cm_npt = [1 1 1; cm_npt];


ax3 = subplot(2,2,3); % total significant patients per square
pcolor(Mptsig'); 
axis equal; axis off; 
% colormap(ax3,raspberry)
% figure('color','w')
cm_sig = cbrewer2('PuRd',7,'seq');
cm_sig = [1 1 1; cm_sig];
colormap(ax3,cm_sig)
% colorbar



ax4 = subplot(2,2,4);  %percent sigificant patients per square
pcolor(Mpercent); 
axis equal; axis off;
cm_percent = cbrewer2('Reds','seq'); %color map
colormap(ax4,cm_percent)

% colormap(ax4,pink_lemonade)

ax2 = subplot(2,2,2); %total patients per square
pcolor(Mnpt'); 
axis equal; axis off; 
% colormap(ax2,mintyfresh)
cm_npt = cbrewer2('Purples',7,'seq');
cm_npt = [1 1 1; cm_npt];
colormap(ax2,cm_npt)




% colormap(spring)




% grid on;
k = 1;

% for l = 1:length(linspace(0,50,7))
