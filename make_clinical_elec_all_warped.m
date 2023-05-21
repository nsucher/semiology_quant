function [emCM,idx_vec]=make_clinical_elec_all_warped(pt)

cd('/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/')
[emC]=getelecs(pt,1);
[emT]=getelecs(pt,2);
[emTM]=getelecs(pt,3);

emCM=nan(size(emC));

nchC=size(emC,1);
nchT=size(emT,1);

idx_vec = zeros(size(emT,1),2); %collect corresponding TDT on 1st col and clinical on 2nd (is same as "num")
idx_vec(:,1) = 1:size(emT);


for i=1:nchC
    dst=[];
    for j=1:nchT
        dst(j)=distance3D(emC(i,:),emT(j,:)); %get distance between clinical electrode i and everryTDT electrode
    end
    mn = min(dst); %find minimum distance, which means that that TDT elkec is the same as the clinical elec of interest
    idx = find(dst==mn); %get the emT row number of it
    if idx > 0
        idx_vec(idx,2) = i;
        emCM(i,:)=emTM(idx,:); %make that emM entry (since each row of emT is the same as the rows of the warped emM version) into the emCM entry (clinical warped version)
    else
        continue
    end
end

d = 1;

