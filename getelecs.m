function [elecmatrix,eleclabels,anatomy]=getelecs(pt,clin1TDT2MNI3)

mainpath='/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/';      %path for OPSCEA folders
%     if exist(mainpath,'dir')==0; mainpath='/Volumes/AN/'; end
%     if ~exist(mainpath,'dir')==7; msgbox('KLEEN_DRIVE not found, you may need to reconnect?'); return; end

ptdir=strcat([mainpath 'OPSCEADATA/'],pt); 
filepath='/Imaging/elecs/';
if clin1TDT2MNI3==1 
    typ='clinical'; 
elseif clin1TDT2MNI3==2
    typ='TDT';
elseif clin1TDT2MNI3==3 
    typ='TDT'; 
end
clinfilename_alternate='clinical_TDT'; 
ext='_elecs_all.mat'; 

cd(mainpath)

if clin1TDT2MNI3<3
    fn=[ptdir filepath typ ext]; 
% elseif exist([ptdir filepath typ '_elecs_all_warped.mat']) == 2
%     fn=[ptdir filepath typ '_elecs_all_warped.mat']; 
else
%     fn=[ptdir filepath typ '_elecs_all2.mat']; 
    fn=[ptdir filepath typ '_elecs_all_warped.mat']; 
end

% elseif exist([ptdir filepath clinfilename_alternate ext]) == 2
fn_alt=[ptdir filepath clinfilename_alternate ext];
% else
%      [emCM,idx_vec]= make_clinical_elec_all_warped(pt);
%      disp('!!')
% end

% fn_ViSAT_alt=['/Volumes/KLEEN_DRIVE/Edwina/ViSAT/ViSAT Data/' pt filepath typ ext];
% fn_ViSAT_alt2=['/Volumes/Edwina/ViSAT/ViSAT Data/' pt filepath typ ext];

if exist(fn)==2;     
    load(fn); 
elseif exist(fn_alt)==2; 
    load(fn_alt); 
% elseif  exist(fn_ViSAT_alt)==2; load(fn_ViSAT_alt); 
% elseif  exist(fn_ViSAT_alt2)==2; load(fn_ViSAT_alt2); 
% elseif exist(elecmatrix) == 2
%     disp(elecmatrix)
else
    [elecmatrix,eleclabels,anatomy]=getOPSCEAelecs(pt,clin1TDT2MNI3); %check OPSCEA pts
    if ~exist('elecmatrix','var'); 
        elecmatrix=[]; eleclabels=[]; anatomy=[]; disp('No result, check path and filename')
    end
end

if ~exist('eleclabels','var'); eleclabels=[]; end
if ~exist('anatomy','var'); anatomy=eleclabels; anatomy=[anatomy cell(size(anatomy,1),1)]; end
