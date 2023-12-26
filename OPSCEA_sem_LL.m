function OPSCEA_sem_LL(uber_ptsz,uber_lat,showlabels,sx_plot,plot_start,plot_end,opscea_path,data_path)
% EXAMPLE USAGE: OPSCEA('UCSF1','01',1,0)

% pt is a string such as 'UCSF4' or 'JaneDoe', acts as a prefix for files below
% sz is a string for '01' or other number of seizure for your patient, acts
% as a secondary prefix (example above becomes UCSF4_01)
% showlabels is:
    % 1 if you want to show the channel labels (default)
    % 0 to hide them AND randomize channels (blinding the reader to the
    % electrode locations of the trace-based ICEEG as in Kleen et al. 2021)
% jumpto allows video to start at a later point, __ seconds ahead

%     Omni-planar and surface casting of epileptiform activity (OPSCEA) (UC
%     Case Number SF2020-281) jointly created by Dr. Jon Kleen, Ben
%     Speidel, Dr. Robert Knowlton, and Dr. Edward Chang is licensed for
%     non-commercial research use at no cost by the Regents of the
%     University of California under CC BY-NC-SA 4.0
%     (https://creativecommons.org/licenses/by-nc-sa/4.0/). Please contact
%     innovation@ucsf.edu if you are interested in using OPSCEA for
%     commercial purposes.

%     The following copyright notice and citation is to be included in any
%     publication, material or media wherein all or a part of Licensed
%     Material is contained, “Certain materials incorporated herein are
%     Copyright © 2016 The Regents of the University of California
%     (REGENTS). All Rights Reserved.

%     Please cite the following paper in your publications if you have used
%     our software in your research, as well as any relevant toolboxes used
%     herein as appropriate (img_pipe, FreeSurfer): Kleen JK, Speidel B,
%     Baud MO, Rao VR, Ammanuel SG, Hamilton LS, Chang EF, Knowlton RC.
%     Accuracy of omni-planar and surface casting of epileptiform activity
%     for intracranial seizure localization. In press at Epilepsia.”

%     Updated by Natalia Sucher Dec 24 2023

if sx_plot
    if ~exist('showlabels','var')||isempty(showlabels);
        showlabels=true; 
    end %default displays ICEEG and depth labels

    % if ~exist('jumpto','var')||isempty(sx_plot); 
    %     sx_plot=0; 
    % end 
    
 %   opscea_path=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
 %   data_path=[opscea_path 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
    if ~exist(data_path,'dir') 
        error('Directory for your data needs to be corrected'); 
    end
    cd(opscea_path);
    
    % ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
    % ptsz = uber_ptsz;
    ptsz = split(uber_ptsz,'_');
    pt = ptsz{1};
    sz = ptsz{2};

    ptpath=[data_path pt '/']; % patient's folder
    szpath= [ptpath uber_ptsz '/']; % specific seizure's folder
    disp(['Running ' pt ', seizure ' sz '...']);
    
    %% Initiate global variables
      global S; % holds general parameters
     % for speed, these are filled during first frame of surfslice then re-used
      global loaf; 
      global sliceinfo; 
      % global I;
    
    
    %% Import parameters
    % for specific seizure 
    [~,prm_allPtSz]=xlsread([data_path 'OPSCEAparams'],'params'); 
        fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
        prm=prm_allPtSz(strcmp(pt,prm_allPtSz(:,1))&strcmp(sz,prm_allPtSz(:,2)),:);
        if isempty(prm) 
            error(['ATTENTION: No entry exists for ' pt ' seizure ' sz ' in the params master sheet']); 
        end
    % Import parameters for patient's specific plot (layout of video frame)
    [~,plt]=xlsread([data_path 'OPSCEAparams'],pt); 
        fields_PLOT=plt(1,:); 
        plt(1,:)=[]; % header for columns of plotting parameters
        plottype=plt(:,strcmpi(fields_PLOT,'plottype')); %type of plot for each subplot (accepts: iceeg, surface, depth, or colorbar)
        %LATERALITY
        % for p = 1:length(plt)
        %     if strcmpi(plottype{p},'cortex'()
        % end
        % plt = plt(1:7,:); % contralateral cortex only
    cd 
    %% prepare subplot specifications
        subplotrow=str2double(plt(:,strcmpi(fields_PLOT,'subplotrow')));
        subplotcolumn=str2double(plt(:,strcmpi(fields_PLOT,'subplotcolumn')));
        subplotstart=plt(:,strcmpi(fields_PLOT,'subplotstart')); 
        subplotstop=plt(:,strcmpi(fields_PLOT,'subplotstop')); 
        for j=1:length(plottype) 
            subplotnum{j,1}=str2double(subplotstart{j}):str2double(subplotstop{j});
        end
        surfaces=plt(:,strcmpi(fields_PLOT,'surfaces'));
        surfacesopacity=plt(:,strcmpi(fields_PLOT,'surfacesopacity'));
        viewangle=lower(plt(:,strcmpi(fields_PLOT,'viewangle')));
    
    %% parcel all individual depth labels, contact #s, and colors. If no depths, make it  =[];
      % depthlabels=plt(:,strcmpi(fields_PLOT,'depthlabels'));
      isdepth=strcmpi(plottype,'depth'); 
      depths=cell(size(isdepth));
      if any(isdepth)
      %   depthEfirst=plt(:,strcmpi(fields_PLOT,'depthEfirst')); 
      %   depthElast=plt(:,strcmpi(fields_PLOT,'depthElast')); 
      %   for j=1:length(depths) 
      %       depths{j}=str2double(depthEfirst{j}):str2double(depthElast{j});
      %   end
      % 
      %   depthcolor=plt(:,strcmpi(fields_PLOT,'depthcolor')); 
      % 
      %   for j=1:length(depthcolor) 
      %       splt=regexp(depthcolor{j},',','split'); 
      %       depthcolor{j}=str2double(splt); 
      %   end 
        pltzoom=str2double(plt(:,strcmpi(fields_PLOT,'pltzoom')));
        pltshowplanes=str2double(plt(:,strcmpi(fields_PLOT,'showplanes')))==1; %logical index of plots in which to show slice planes
      end
    
    %% Get time segments within the ICEEG file to use
        VIDstart=prm(:,strcmpi(fields_SZ,'VIDstart')); VIDstop=prm(:,strcmpi(fields_SZ,'VIDstop')); %chunk of data (seconds into ICEEG data file) to use from the whole ICEEG data clip for the video
        S.VIDperiod=[str2double(VIDstart{1}) str2double(VIDstop{1})];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vid_period = S.VIDperiod; %NS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BLstart=prm(:,strcmpi(fields_SZ,'BLstart')); BLstop=prm(:,strcmpi(fields_SZ,'BLstop')); %chunk of data (seconds into ICEEG data file) to use for baseline (for z-score step)
        S.BLperiod=[str2double(BLstart{1}) str2double(BLstop{1})];
    
    %transform, scaling, and display options
    S.llw=str2double(prm{strcmp('llw',fields_SZ)}); %default linelength window (in seconds)
    S.iceeg_scale=prm{strcmp('iceeg_scale',fields_SZ)}; %percentile (number >50 and <100), used here similar to gain ICEEG waveform display, usually 95
        if ischar(S.iceeg_scale); S.iceeg_scale=str2double(S.iceeg_scale); end 
    S.fps=str2double(prm{strcmp('fps',fields_SZ)});             %frames per sec of ICEEG (default 15)
    S.cax=str2double(regexp(prm{strcmp('cax',fields_SZ)},',','split'));         %color axis for heatmap
    S.gsp=str2double(prm{strcmp('gsp',fields_SZ)}); %gaussian spreading parameter (default 10)
        params={'iceeg_scale','fps','cax','gsp'}; 
        paramsnans=isnan([(isnan(S.iceeg_scale) | S.iceeg_scale<=50 | S.iceeg_scale>=100)   S.fps   any(isnan(S.cax)) S.gsp]); 
        if any(paramsnans); error(['ATTENTION OPSCEA USER: The "' params{paramsnans} '" term(s) is/are in an incorrect format (perhaps number instead of string), check excel seizure parameter sheet']); 
        end 
      cm=prm{strcmp('cm',fields_SZ)};
      switch cm; case 'cmOPSCEAcool'; cm=cmOPSCEAcool; 
                 case 'cmOPSCEAjet'; cm=cmOPSCEAjet; 
      end
    S.cm=cm; %colormap to use for heatmap
    S.iceegwin=str2double(prm{strcmp('iceegwin',fields_SZ)}); %how much trace-based ICEEG to view at a time in the ICEEG window
    S.marg=str2double(prm{strcmp('marg',fields_SZ)}); %offset of real-time LL txform from beginning of viewing window (in sec; converts to samples below)
    S.slicebright=str2double(prm{strcmp('slicebright',fields_SZ)}); if isnan(S.slicebright); S.slicebright=0; end %brighten up slices (usually 0 to 50)
    
    
    % additional adjustment for display window
    S.VIDperiod=[S.VIDperiod(1)-S.marg   S.VIDperiod(2)+S.iceegwin-S.marg]; 
    S.fields=fields_SZ; clear fields
    
    S.prm=prm; clear prm
    S.prm_allPtSz=prm_allPtSz; clear prm_allPtSz
    S=orderfields(S); %alphabetize the structure fields for ease of use/search
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S.sliceplane='c'; % calculate omni-planar slice angles with respect to coronal (c) plane
    
    %% load ICEEG data, and the bad channels verified for that specific data
    load([szpath uber_ptsz])
    load([szpath uber_ptsz '_badch']); 
    % rename and clear old format of electrode files -NS
    if exist('ppEEG','var')
        d = ppEEG; clear ppEEG; end;
    if exist('fs','var')
        sfx = fs; clear fs; end;
    if exist('bad_chs','var')
        badch = bad_chs; clear bad_chs; end; 
    % end of edit -NS
    if size(d,1)>size(d,2); d=d'; end % orient to channels by samples
    [~,ntp]=size(d); f=1;  %EDITED JK 5/4
    % disp(['Length of data to play for video is ' num2str(round(ntp/sfx)) 'sec'])
    
    % error checks for selected time periods
    if any([S.VIDperiod(1) S.BLperiod(1)]<0)
        error('VIDperiod is out of bounds of file (time < 0). Check both VIDstart and BLstart times and make sure the "marg" value (subtracted from VIDstart and BLstart), cannot be < 0'); 
    elseif any([S.VIDperiod(2) S.BLperiod(2)]>ntp)
        error('VIDperiod is beyond the length of the file. Check VIDstop and BLstop times vs. length of actual ICEEG data'); 
    end 
    
    %% locate and load electrode file for labels and XYZ coordinates
    if exist([ptpath 'Imaging/Elecs/Electrodefile.mat'])
        load([ptpath 'Imaging/Elecs/Electrodefile.mat']); 
    elseif exist([ptpath 'Imaging/elecs/clinical_elecs_all.mat']) % access variables in old format NS
         load([ptpath 'Imaging/elecs/clinical_elecs_all.mat']); % NS
    end
    if ~exist('anatomy','var') 
        anatomy=cell(size(elecmatrix,1),4); 
    end
    if size(anatomy,1)>size(elecmatrix,1) 
        anatomy(size(elecmatrix,1)+1:end,:)=[]; 
    end
    anat=anatomy; 
    clear anatomy; 
    if size(anat,2)>size(anat,1)
        anat=anat'; 
    end
    if size(anat,2)==1
        anat(:,2)=anat(:,1); 
    end

    if ~exist('eleclabels','var')
        eleclabels=anat(:,1); 
    end
   em=elecmatrix; 
   clear elecmatrix; 
   emnan=isnan(mean(em,2)); 
   em(emnan,:)=0; 
   EKGorREF=strcmpi('EKG1',anat(:,1))|strcmpi('EKG2',anat(:,1))|strcmpi('EKG',anat(:,2))|strcmpi('EKGL',anat(:,2))|strcmpi('REF',anat(:,1)); anat(EKGorREF,:)=[]; em(EKGorREF,:)=[]; eleclabels(EKGorREF,:)=[]; 
   d(size(anat,1)+1:end,:)=[];
   badch(size(anat,1)+1:end)=[];
   [nch,~]=size(d); %%%%%%%% edited JK 5/4
    
    isR=nansum(em(:,1))>0; isL=isR~=1; %handy binary indicators for laterality
    
    %% load meshes you want to plot
    meshpath='Imaging/Meshes/';
    Rcortex=load([ptpath meshpath pt '_rh_pial.mat']); loaf.rpial=Rcortex; Rcrtx=Rcortex.cortex; clear Rcortex
    Lcortex=load([ptpath meshpath pt '_lh_pial.mat']); loaf.lpial=Lcortex; Lcrtx=Lcortex.cortex; clear Lcortex
    
    for i=1:length(surfaces); hippentry(i)=~isempty(strfind(surfaces{i},'hipp')); amygentry(i)=~isempty(strfind(surfaces{i},'amyg')); end; 
    errmsg='ATTN: MISSING A MESH, need to add this mesh file to directory (or remove/omit from frame): ';
      if any(hippentry); Rhipp=[ptpath meshpath 'subcortical/rHipp_subcort.mat']; Lhipp=Rhipp; Lhipp(end-16)='l'; if exist(Rhipp,'file'); Rhipp=load(Rhipp); Rhipp=Rhipp.cortex;  Lhipp=load(Lhipp); Lhipp=Lhipp.cortex;     else; error([errmsg 'hipp']); end; end
      if any(amygentry); Ramyg=[ptpath meshpath 'subcortical/rAmgd_subcort.mat']; Lamyg=Ramyg; Lamyg(end-16)='l'; if exist(Ramyg,'file'); Ramyg=load(Ramyg); Ramyg=Ramyg.cortex;  Lamyg=load(Lamyg); Lamyg=Lamyg.cortex;     else; error([errmsg 'amyg']); end; end
    drows=find(strcmp(plottype,'depth'))'; ndepths=length(drows);
    
    % depthch=[]; 
    % for i=1:length(drows); 
    %     depthch=[depthch depths{drows(i)}]; 
    % end; clear i %identify all depth electrode channels
    
    %% get xyz limits for plotting purposes
    perim=1; % how many millimeters away from brain/electrodes boundaries to set the colorcoded plane perimeter, recommend >0 to avoid skimming brain surface (default 1mm)
    axl(:,:,1)=[min([Rcrtx.vert; Lcrtx.vert]); max([Rcrtx.vert; Lcrtx.vert])]; %min and max of MESH VERTICES' x y z coordinates (2x3)
    axl(:,:,2)=[min(em); max(em)]; %%min and max of ELECTRODES' x y z coordinates (2x3) and stack them (2x3x2)
    axl=[min(axl(1,:,:),[],3); max(axl(2,:,:),[],3)]; % Get the minima and maxima of both
    axislim=reshape(axl,1,6)+[-1 1 -1 1 -1 1]*perim; clear axl %Use the, to define the axis boundaries, and add additional perimeter (perim)
    
    %% formatting checks, and consolidation of bad channels
    ns=unique( [find(badch);   find(isnan(mean(em,2)));   find(isnan(mean(d,2)))]  ); % bad channels: those that are pre-marked, or if NaNs in their coordinates or ICEEG data traces
    nns=true(nch,1); 
    nns(ns)=0; %nns=find(nns); %consolidate bad channels and those with NaNs
    
    %remove data channels previously clipped at end. Only include that which has electrode coordinates (intracranial)
    if size(em,1)>size(d,1); nch=size(em,1); d(nch+1:end,:)=[]; LL(nch+1:end,:)=[]; nns(nch+1:end)=[]; ns(ns>size(em,1))=[]; 
       fprintf(2, 'ALERT: Clipping off extra bad channel entries (make sure you have the right ICEEG and bad channel files loaded)\n');
    end
    
    
    %% ICEEG data processing and transform
    
    % Filter out < 1 Hz (and up to nyquist) out to help decrease movement
    % artifact and improve stability during ICEEG trace plotting
    [b,a]=butter(2,[1 round(sfx/2-1)]/(sfx/2),'bandpass'); 
    d(nns,:) = filtfilt(b,a,d(nns,:)')';
    
    % Line-length transform
    L=round(S.llw*sfx)-1; % number of samples to calculate line length
    LL=nan(size(d)); 
    for i=1:size(d,2)-L; LL(:,i)=sum(abs(diff(d(:,i:i+L),1,2)),2); end
    
    %Normalize LL (to baseline period) as z-scores, "zLL"
    BLstartsample=max([sfx*S.BLperiod(1) 1]); BLendsample=sfx*S.BLperiod(2); 
    for i=1:nch % z-score using channel-specific baseline
        LL(i,:)=(LL(i,:)-nanmean(LL(i,BLstartsample:BLendsample)))/nanstd(LL(i,BLstartsample:BLendsample)); 
    end 
    
    % make a timtestamp vector
    ts=0:1/sfx:size(d,2)*(1/sfx)-1/sfx;
    
    
    % Extract the period of data to be used for the video (remove flanking data)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if S.VIDperiod(2)*sfx > size(d,2)
        vidperiodidx=round(S.VIDperiod(1)*sfx+1):size(d,2);
    else
        vidperiodidx=round(S.VIDperiod(1)*sfx+1):S.VIDperiod(2)*sfx;
    end

    d=d(:,vidperiodidx);
    ntp=length(vidperiodidx);
    LL=LL(:,vidperiodidx); 
    ts=ts(vidperiodidx); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Scaling ICEEG and setting windows for simultaneous trace-based display
    S.iceeg_scale=S.iceeg_scale+(100-S.iceeg_scale)/2; S.iceeg_scale=[100-S.iceeg_scale S.iceeg_scale]; %conversion to two-tailed percentile
        scl=2/diff(prctile(reshape(d,1,numel(d)),S.iceeg_scale)); 
    S.marg=round(S.marg*sfx); %offset of real-time LL txform from beginning of viewing window
%     jumpto=S.marg+round(jumpto*sfx); %Jump ahead (sec), so video starts this much farther into the file if desired. Input argument.
    sx_plot=S.marg+round(sx_plot*sfx)-sfx; %Jump ahead (sec), so video starts this much farther into the file if desired. Input argument.

    S.fram=round(sfx/S.fps);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLOTTING TIME!
    sliceinfo=[]; loaf.vrf=[]; loaf.apasrf=[]; loaf.normloaf=[]; sliceinfo.viewangle=zeros(size(plt,1),3); sliceinfo.azel=[]; sliceinfo.corners=[]; loaf.isR=isR; loaf.isL=isL; 
    clear F; 
    ytl=eleclabels(nns,1); 
    nch=length(find(nns)); 
    
    chanorder=1:size(d(nns,:),1); if ~showlabels; chanorder=randperm(size(d(nns,:),1)); end % if desired, blinds user by randomizing channel order
    figure('color','w','Position',[1 5 1280 700]); 
    % frametimpoints=jumpto:S.fram:ntp-sfx*S.iceegwin; % timepoint index of each frame to be rendered
    round_sx_plot = floor(sx_plot);
    for i = round_sx_plot
        subplot(1,1,1); %clears all axes, to start fresh each frame
        w8s=LL(:,i); 
        w8s(~nns)=0; %make weights for electrodes, and set NaNs (bad channels) to zero
        for j=1:size(plt,1)
            subplot(subplotrow(j),subplotcolumn(j),subplotnum{j}); %to do: MAKE ONLY ONE BRAIN/SUBPLOT (keeping zLL plot and iceeg plot)
        switch upper(plottype{j,1})
          case 'ECOG' %plot the raw data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
               %subplot(1,400,[330:350]) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dtoplot=d(nns,((i-S.marg+1):(i-S.marg+1)+sfx*S.iceegwin))-1;
                tstoplot=ts((i-S.marg+1):(i-S.marg+1)+sfx*S.iceegwin)-1;
                % dtoplot=d(nns,(i-S.marg):(i-S.marg)+sfx*S.iceegwin);
                % tstoplot=ts((i-S.marg):(i-S.marg)+sfx*S.iceegwin);

                shift = repmat(-1*(1:nch)',1,size(dtoplot,2)); 
                plot(tstoplot,dtoplot*scl+shift,'k','LineWidth',.000001);
                  ylim([-nch-1 0])
                  axis tight; xlabel('time (sec)')
                  hold on;
                fill((ones(1,4)*ts(i)+[0 S.llw S.llw 0])-1,[.5 .5 -nch-1.25 -nch-1.25],[.75 0 .75],'facealpha',.25,'edgealpha',1); hold off; % overlay transform window
                  xlabel('Time (seconds)'); 
                  textfactor=min([ceil((length(find(nns))-80)/20) 4]); %scale text size
                  if showlabels; set(gca,'ytick',-length(ytl):-1,'yticklabel',flipud(ytl(chanorder)),'fontsize',8-textfactor); else; set(gca,'ytick',[]); ylabel('Channels (randomized order)'); end
                  set(gca,'ylim',[-(nch)-1.25 .25])
%                   text(i/sfx-.01+S.llw*.36,2,'v'); text(repmat(i/sfx+.01705+S.llw*.4,1,4),2.5:1:5.5,{'|','|','|','|'}) %draws an arrow pointing to the transform window
                  ttl1=title('ICEEG'); set(ttl1,'fontsize',10)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
              hold on;           
              % plot([sx_plot sx_plot],ylim,'m-') % vertical magenta line that follows time
              % hold on;
            %Clear non-labeled channels beyond size of number of electrode rows
            %Vectorize unwanted channels 
              noneed=false(size(anat,1),4);
              anat_rows = size(anat,1); 
              badch = badch(1:anat_rows); % cut down badch to eliminate extra bad channels after anatomy rows size
              for num_rows=1:size(anat,1)
                  noneed(num_rows,1) = contains(lower(anat{num_rows,4}),'ctx'); % cell array containing row index of strings with "ctx" in u1
                  noneed(num_rows,2) = contains(lower(anat{num_rows,4}),'wm'); % cell array containing row index of strings with "wm" in u1
                  noneed(num_rows,3) = contains(lower(anat{num_rows,4}),'white-matter'); % cell array containing row index of strings with "Right-Cerebral-White-Matter" in u1             
                  noneed(num_rows,4) = contains(lower(anat{num_rows,4}),'unknown'); % cell array containing row index of strings with "Unknown" in u1
                  noneed(num_rows,5) = contains(lower(anat{num_rows,4}),'vent'); % cell array containing row index of strings with "Unknown" in u1             
              end
              noneed=any(noneed,2);
    
              noneed = noneed | badch; % now is either uncessary (noneed) or bad channels
    

              new_LL=LL; %do not clear LL, used to construct w8s
              new_LL(noneed,:)=[];

              new_anat=anat; %do not clear anat, used to construct noneed 
              new_anat(noneed,:)=[];

              new_em=em; %do not clear anat, used to construct noneed 
              new_em(noneed,:)=[];


              new_w8s = w8s;
              new_w8s(noneed,:) = [];
              subplot(2,100,62:100)
              LL_plot(new_anat,new_LL,ts,sx_plot,plot_start,plot_end,sfx,S.cax);

             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
          case 'COLORBAR' 
                    colormap_ns(gca,S.cm(ceil(size(S.cm,1)/2):end,:)); %Z-scores above baseline

                    hold on; plot(1,1); title(''); axis off; 

                    cb=colorbar; 

                    cb.Ticks=[.5 1]; cb.Limits=[.5 1];
                    cb.TickLabels={0,num2str(S.cax(2))};
                    cb.FontSize=11; cb.Location='west'; 
% ylabel(cb,'z-scores','fontsize',12);


          case 'SURFACE' % plotting surfaces only
              if matches(surfaces{j},'rcortex')
                  hold off; 
                  srf=regexp(surfaces{j},',','split'); % list the specific surfaces wanted for this subplot
                  srfalpha=regexp(surfacesopacity{j},',','split'); % list their corresponding opacities (values from 0 to 1; 0=invisible, 1=opaque)
                  if length(srf)~=length(srfalpha)
                      msgbox('Number of surface to plot does not match number of alpha designations, check excel sheet'); 
                      return; 
                  end
                  acceptedterms={'rcortex','lcortex','rhipp','lhipp','ramyg','lamyg','wholebrain'};
                    for s=1:length(srf)
                        srf{s}=lower(srf{s}); %convert to lower case for easier string matching
                      if ~isempty(intersect(srf{s},acceptedterms)) %make sure user specified accepted terminologies for the meshes
                        switch srf{s}; %see below for case "wholebrain"s
                            case 'rcortex'; srfplot=Rcrtx; 
                            case 'lcortex'; srfplot=Lcrtx; 
                            case 'rhipp';   srfplot=Rhipp; 
                            case 'lhipp';   srfplot=Lhipp; 
                            case 'ramyg';   srfplot=Ramyg; 
                            case 'lamyg';   srfplot=Lamyg; 
                        end
                      end

                      % plot the individual heatmapped surface
                      cd(opscea_path)
                      if exist('srfplot','var') 
                          hh=ctmr_gauss_plot_edited(srfplot,new_em,new_w8s,S.cax,0,S.cm,S.gsp); 
                          alpha(hh,str2double(srfalpha{s})); % Adjust opacity specified for that row
                      else 
                          disp(['ALERT: One of the entries in row ' num2str(j) ' is not a valid entry, accepts:']); 
                          disp(acceptedterms); 
                      end
                    end
                    if isempty(intersect(srf{s},{'rcortex','lcortex'}))||strcmpi(srf,'wholebrain') %for glass brain (hipp and/or amyg only) and wholebrain plots
                        glass1=ctmr_gauss_plot_edited(Rcrtx,new_em,new_w8s,S.cax,0,S.cm,S.gsp); alpha(glass1,.1); 
                        glass2=ctmr_gauss_plot_edited(Lcrtx,new_em,new_w8s,S.cax,0,S.cm,S.gsp); alpha(glass2,.1); 
                        if ~pltshowplanes(j) 
                            plot3(new_em(:,1),new_em(:,2),new_em(:,3),'o','Color','k','MarkerFaceColor','k','markersize',2.5); 
                        end 
                    else 
                        plot3(new_em(:,1),new_em(:,2),new_em(:,3),'o','Color','k','MarkerFaceColor','k','markersize',2.5); %plot electrodes
                    end
                    cameratoolbar('setmode',''); 
                    litebrain(viewangle{j},.9); 
                    wb=strcmpi(srf,'wholebrain'); 
                    if any(wb) 
                        alpha(glass1,srfalpha{wb}); 
                        alpha(glass2,srfalpha{wb}); 
                    end
                    if strcmp(viewangle{j},'i')
                        view(90+isL*180,270)
                    end
                    if pltshowplanes(j)||(strcmp(viewangle{j},'i')&&~isempty(intersect(srf,'wholebrain'))) 
                        view(180,270); 
                    end %orients the "show planes" slice to a classic axial perspective
                    axis(axislim); 
                    if strcmpi(viewangle{j},'i')||strcmpi(viewangle{j},'s')||strcmpi(viewangle{j},'a')||strcmpi(viewangle{j},'p');  
                        if ~strcmpi(pt,'NO181') 
                            axis([axislim(1)*isL+10*isR axislim(2)*isR+10*isL  axislim(3:6)]); 
                        end 
                    end
                    zoom(pltzoom(j)); 
                    hold on;                         
                    cd(opscea_path)
                    colormap_ns(gca,S.cm); set(gca,'Clipping','off')
                    clear srfplot
              end
            end
        end  
     end
   end
end
