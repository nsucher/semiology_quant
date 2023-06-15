function [elecmatrix,eleclabels,anatomy]=getOPSCEAelecs(pt,clin1TDT2MNI3)

mainpath='/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/';      %path for OPSCEA folders
%     if exist(mainpath,'dir')==0; mainpath='/Volumes/OPSCEA/'; end
%     if ~exist(mainpath,'dir')==7; msgbox('KLEEN_DRIVE not found, you may need to reconnect?'); return; end

ptdir=strcat([mainpath 'OPSCEADATA/'],pt); 
filepath='/Imaging/elecs/';
if clin1TDT2MNI3==1; typ='clinical'; elseif clin1TDT2MNI3==2; typ='TDT'; elseif clin1TDT2MNI3==3; typ='TDT'; end;
clinfilename_alternate='clinical_TDT'; 
ext='_elecs_all.mat'; 


if clin1TDT2MNI3<3; fn=[ptdir filepath typ ext]; else; fn=[ptdir filepath typ '_elecs_all_warped.mat']; end
fn_alt=[ptdir filepath clinfilename_alternate ext];

if      exist(fn)==2;     
    load(fn); 
elseif  exist(fn_alt)==2; load(fn_alt); 
else elecmatrix=[]; eleclabels=[]; anatomy=[]; disp('No result, check path and filename')
end

if ~exist('eleclabels','var'); eleclabels=[]; end
if ~exist('anatomy','var'); anatomy=eleclabels; anatomy=[anatomy cell(size(anatomy,1),1)]; end
