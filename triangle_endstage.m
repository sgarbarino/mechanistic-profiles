
indAD = [nodAD_ES';troAD_ES';traAD_ES';proAD_ES'];
indMS = [nodMS_ES';troMS_ES';traMS_ES';proMS_ES'];
indRS = [nodRS_ES';troRS_ES';traRS_ES';proRS_ES'];
% distances
di_ADtoAD = pdist2(indAD',centerAD_ES);
di_ADtoMS = pdist2(indAD',centerMS_ES);
di_ADtoRS = pdist2(indAD',centerRS_ES);

di_MStoAD = pdist2(indMS',centerAD_ES);
di_MStoMS = pdist2(indMS',centerMS_ES);
di_MStoRS = pdist2(indMS',centerRS_ES);

di_RStoAD = pdist2(indRS',centerAD_ES);
di_RStoMS = pdist2(indRS',centerMS_ES);
di_RStoRS = pdist2(indRS',centerRS_ES);

diADtocenters = [di_ADtoAD di_ADtoMS di_ADtoRS];
diMStocenters = [di_MStoAD di_MStoMS di_MStoRS];
diRStocenters = [di_RStoAD di_RStoMS di_RStoRS];

% from distances to weights

weiADtocenters = 1./diADtocenters;
weiADtocenters(:,1) = -.5+ weiADtocenters(:,1);
weiADtocenters(1:400,2) =  abs(weiADtocenters(1:400,2) -0); %%
weiADtocenters(1:400,3) =  abs(weiADtocenters(1:400,3) -1);
weiADtocenters(401:800,2) =  abs(weiADtocenters(401:800,2) -1);
weiADtocenters(401:800,3) =  abs(weiADtocenters(401:800,3) -0);
weiADtocenters(801:end,2) =  abs(weiADtocenters(801:end,2) -0);
weiADtocenters(801:end,3) =  abs(weiADtocenters(801:end,3) -0);
weiADtocenters = weiADtocenters./sum(weiADtocenters,2);

weiMStocenters = 1./diMStocenters;
weiMStocenters(:,2) = -.5+ weiMStocenters(:,2);
weiMStocenters(1:14,1) =  abs(weiMStocenters(1:14,1) -1);
weiMStocenters(1:14,3) =  abs(weiMStocenters(1:14,3) -0);
weiMStocenters(15:30,1) =  abs(weiMStocenters(15:30,1) -0.5); %%
weiMStocenters(15:30,3) =  abs(weiMStocenters(15:30,3) -1);
weiMStocenters(31:end,1) =  abs(weiMStocenters(31:end,1) -0);
weiMStocenters(31:end,3) =  abs(weiMStocenters(31:end,3) -1);

weiMStocenters = weiMStocenters./sum(weiMStocenters,2);

weiRStocenters = 1./diRStocenters;
weiRStocenters(1:1000,2) =  abs(weiRStocenters(1:1000,2) -1);
weiRStocenters(1:1000,1) =  abs(weiRStocenters(1:1000,1) -0);
weiRStocenters(1001:2000,2) =  abs(weiRStocenters(1001:2000,2) -0);
weiRStocenters(1001:2000,1) =  abs(weiRStocenters(1001:2000,1) -1);
weiRStocenters(2001:end,2) =  abs(weiRStocenters(2001:end,2) -1);
weiRStocenters(2001:end,1) =  abs(weiRStocenters(2001:end,1) -1);

weiRStocenters = weiRStocenters./sum(weiRStocenters,2);

%%

[h,~,~]=terplot;
[hter_RS,xRS,yRS]=ternaryc(weiRStocenters(:,1),weiRStocenters(:,2),weiRStocenters(:,3));
[hter_AD,xAD,yAD]=ternaryc(weiADtocenters(:,1),weiADtocenters(:,2),weiADtocenters(:,3));
[hter_MS,xMS,yMS]=ternaryc(weiMStocenters(:,1),weiMStocenters(:,2),weiMStocenters(:,3));
set(hter_AD,'marker','o','color','k','markerfacecolor','r','markersize',8)
set(hter_MS,'marker','o','color','k','markerfacecolor','g','markersize',8)
set(hter_RS,'marker','o','color','k','markerfacecolor','b','markersize',8)


dADMS = pdist2([xAD;yAD]', [0;1]');
dADRS = pdist2([xAD;yAD]', [0.5;sqrt(1-(0.5)^2)]');
dADAD = pdist2([xAD;yAD]', [0;0]');

dMSMS = pdist2([xMS;yMS]', [1;0]');
dMSRS = pdist2([xMS;yMS]', [0.5;sqrt(1-(0.5)^2)]');
dMSAD = pdist2([xMS;yMS]', [0;0]');

dRSMS = pdist2([xRS;yRS]', [1;0]');
dRSRS = pdist2([xRS;yRS]', [0.5;sqrt(1-(0.5)^2)]');
dRSAD = pdist2([xRS;yRS]', [0;0]');

diADC = [dADAD  dADMS  dADRS];
diMSC = [dMSAD  dMSMS  dMSRS];
diRSC = [dRSAD  dRSMS  dRSRS];

[~,tAD] = min(diADC');
idxADout = find(tAD~=1);
[~,tMS] = min(diMSC');
idxMSout = find(tMS~=2);
[~,tRS] = min(diRSC');
idxRSout = find(tRS~=3);


set(hter_AD(idxADout),'marker','d','MarkerEdgeColor','r','linewidth',2,'markerfacecolor','k','markersize',20) 
set(hter_MS(idxMSout),'marker','d','MarkerEdgeColor','g','linewidth',2,'markerfacecolor','k','markersize',20) 
set(hter_RS(idxRSout),'marker','d','MarkerEdgeColor','b','linewidth',2,'markerfacecolor','k','markersize',20) 

hlabels=terlabel('AD','PPMS','HA');
set(hlabels,'Linewidth',30)

%% outliers analysis - just if access to data
% % AD analysis
% [~,tmpAD] = min(diADtocenters');
% idxAD_inMS = find(tmpAD==2);
% idxAD_inRS = find(tmpAD==3);
% 
% AD_inMS = demoFINAL_RID_PAT_AD(idxAD_inMS,:); 
% AD_inRS = demoFINAL_RID_PAT_AD(idxAD_inRS,:); 
% 
% remaining_AD_fromMS = demoFINAL_RID_PAT_AD; 
% remaining_AD_fromMS(idxAD_inMS,:)=[]; 
% remaining_AD_fromRS = demoFINAL_RID_PAT_AD;
% remaining_AD_fromRS(idxAD_inRS,:)=[]; 
% 
% % MS analysis
% [~,tmpMS] = min(diMStocenters');
% idxMS_inAD = find(tmpMS==1);
% idxMS_inRS = find(tmpMS==3);
% 
% MS_inAD = demoFINAL_RID_PPMS(idxMS_inAD,:);
% MS_inRS = demoFINAL_RID_PPMS(idxMS_inRS,:);
% 
% remaining_MS_fromAD = demoFINAL_RID_PPMS;
% remaining_MS_fromAD(idxMS_inAD,:)=[]; 
% remaining_MS_fromRS = demoFINAL_RID_PPMS;
% remaining_MS_fromRS(idxMS_inRS,:)=[]; 
% 
% % RS analysis
% % [~,tmpRS] = min(distRStocenters');
% % idxRS_inAD = find(tmpRS==1);
% % idxRS_inMS = find(tmpRS==2);
% 
% %%%% new data analysis Data_demo 
% load new_RS_data.mat
% % RS_inAD = Data_demo(idxRS_inAD,:);
% % RS_inMS = Data_demo(idxRS_inMS,:);
% % 
% % remaining_RS_fromAD = Data_demo;
% % remaining_RS_fromAD(idxRS_inAD,:)=[]; 
% % remaining_RS_fromMS = Data_demo;
% % remaining_RS_fromMS(idxRS_inMS,:)=[]; 
% RS_inAD = Data_crosssec_764_closerAD_ES;
% RS_inMS = Data_crosssec_1529_closerMS_ES;
% remaining_RS_fromAD = Data_crosssec_all_clean;
% remaining_RS_fromAD(idxRS_inAD_ES,:)=[]; 
% remaining_RS_fromMS = Data_crosssec_all_clean;
% remaining_RS_fromMS(idxRS_inMS_ES,:)=[]; 
% 
% 
% 
% mode(AD_inRS.DX_bl)
% mode(remaining_AD_fromRS.DX_bl)
% % [H,P]=ttest2(AD_inRS.DX_bl,remaining_AD_fromRS.DX_bl)
% mode(AD_inRS.AGE)
% mode(remaining_AD_fromRS.AGE)
% [H,P]=ttest2(AD_inRS.AGE,remaining_AD_fromRS.AGE)
% mode(AD_inRS.PTEDUCAT)
% mode(remaining_AD_fromRS.PTEDUCAT)
% [H,P]=ttest2(AD_inRS.PTEDUCAT,remaining_AD_fromRS.PTEDUCAT)
% mode(AD_inRS.PTGENDER)
% mode(remaining_AD_fromRS.PTGENDER)
% [H,P]=ttest2(AD_inRS.PTGENDER,remaining_AD_fromRS.PTGENDER)
% mean(AD_inRS.MMSE_bl)
% mean(remaining_AD_fromRS.MMSE_bl)
% [H,P]=ttest2(AD_inRS.MMSE_bl,remaining_AD_fromRS.MMSE_bl)
% nanmean(AD_inRS.ADAS13_bl)
% nanmean(remaining_AD_fromRS.ADAS13_bl)
% [H,P]=ttest2(AD_inRS.ADAS13_bl,remaining_AD_fromRS.ADAS13_bl)
% nanmean(AD_inRS.AV45_bl)
% nanmean(remaining_AD_fromRS.AV45_bl)
% [H,P]=ttest2(AD_inRS.AV45_bl,remaining_AD_fromRS.AV45_bl) %%%% differences
% nanmean(AD_inRS.ABETA_UPENNBIOMK9_04_19_17)
% nanmean(remaining_AD_fromRS.ABETA_UPENNBIOMK9_04_19_17)
% [H,P]=ttest2(AD_inRS.ABETA_UPENNBIOMK9_04_19_17,remaining_AD_fromRS.ABETA_UPENNBIOMK9_04_19_17)
% nanmean(AD_inRS.PTAU_UPENNBIOMK9_04_19_17)
% nanmean(remaining_AD_fromRS.PTAU_UPENNBIOMK9_04_19_17)
% [H,P]=ttest2(AD_inRS.PTAU_UPENNBIOMK9_04_19_17,remaining_AD_fromRS.PTAU_UPENNBIOMK9_04_19_17)
% mode(AD_inRS.APOE4)
% mode(remaining_AD_fromRS.APOE4)
% [H,P]=ttest2(AD_inRS.APOE4,remaining_AD_fromRS.APOE4)
% nanmean(AD_inRS.CDRSB_bl)
% nanmean(remaining_AD_fromRS.CDRSB_bl)
% [H,P]=ttest2(AD_inRS.CDRSB_bl,remaining_AD_fromRS.CDRSB_bl)
% %% comparisons AD in MS vs AD not in MS
% 
% mode(AD_inMS.DX_bl)
% mode(remaining_AD_fromMS.DX_bl)
% % [H,P]=ttest2(AD_inRS.DX_bl,remaining_AD_fromRS.DX_bl)
% mode(AD_inMS.AGE)
% mode(remaining_AD_fromMS.AGE)
% [H,P]=ttest2(AD_inMS.AGE,remaining_AD_fromMS.AGE)
% mode(AD_inMS.PTEDUCAT)
% mode(remaining_AD_fromMS.PTEDUCAT)
% [H,P]=ttest2(AD_inMS.PTEDUCAT,remaining_AD_fromMS.PTEDUCAT)
% mode(AD_inMS.PTGENDER)
% mode(remaining_AD_fromMS.PTGENDER)
% [H,P]=ttest2(AD_inMS.PTGENDER,remaining_AD_fromMS.PTGENDER)
% mean(AD_inMS.MMSE_bl)
% mean(remaining_AD_fromMS.MMSE_bl)
% [H,P]=ttest2(AD_inMS.MMSE_bl,remaining_AD_fromMS.MMSE_bl)
% nanmean(AD_inMS.ADAS13_bl)
% nanmean(remaining_AD_fromMS.ADAS13_bl)
% [H,P]=ttest2(AD_inMS.ADAS13_bl,remaining_AD_fromMS.ADAS13_bl)
% nanmean(AD_inMS.AV45_bl)
% nanmean(remaining_AD_fromMS.AV45_bl)
% [H,P]=ttest2(AD_inMS.AV45_bl,remaining_AD_fromMS.AV45_bl) %%%% differences
% nanmean(AD_inMS.ABETA_UPENNBIOMK9_04_19_17)
% nanmean(remaining_AD_fromMS.ABETA_UPENNBIOMK9_04_19_17)
% [H,P]=ttest2(AD_inMS.ABETA_UPENNBIOMK9_04_19_17,remaining_AD_fromMS.ABETA_UPENNBIOMK9_04_19_17)
% nanmean(AD_inMS.PTAU_UPENNBIOMK9_04_19_17)
% nanmean(remaining_AD_fromMS.PTAU_UPENNBIOMK9_04_19_17)
% [H,P]=ttest2(AD_inMS.PTAU_UPENNBIOMK9_04_19_17,remaining_AD_fromMS.PTAU_UPENNBIOMK9_04_19_17)
% mode(AD_inMS.APOE4)
% mode(remaining_AD_fromMS.APOE4)
% [H,P]=ttest2(AD_inMS.APOE4,remaining_AD_fromMS.APOE4)
% nanmean(AD_inMS.CDRSB_bl)
% nanmean(remaining_AD_fromMS.CDRSB_bl)
% [H,P]=ttest2(AD_inMS.CDRSB_bl,remaining_AD_fromMS.CDRSB_bl)
% 
% %% comparison MS in AD vs MS not in AD
% nanmean(cell2mat(MS_inAD.EDSS))
% nanmean(cell2mat(remaining_MS_fromAD.EDSS))
% [H,P]=ttest2(cell2mat(MS_inAD.EDSS),cell2mat(remaining_MS_fromAD.EDSS))
% 
% nanmean(MS_inAD.age_baseline)
% nanmean(remaining_MS_fromAD.age_baseline)
% [H,P]=ttest2(MS_inAD.age_baseline,remaining_MS_fromAD.age_baseline)
% 
% nanmean(cell2mat(MS_inAD.disease_duration))
% nanmean(cell2mat(remaining_MS_fromAD.disease_duration))
% [H,P]=ttest2(cell2mat(MS_inAD.disease_duration),cell2mat(remaining_MS_fromAD.disease_duration))
% 
% mode((MS_inAD.gender))
% mode((remaining_MS_fromAD.gender))
% [H,P]=ttest2((MS_inAD.gender),(remaining_MS_fromAD.gender))
% 
% 
% %% comparisons MS in HA vs MS not in HA
% nanmean(cell2mat(MS_inRS.EDSS))
% nanmean(cell2mat(remaining_MS_fromRS.EDSS))
% [H,P]=ttest2(cell2mat(MS_inRS.EDSS),cell2mat(remaining_MS_fromRS.EDSS))
% 
% nanmean(MS_inRS.age_baseline)
% nanmean(remaining_MS_fromRS.age_baseline)
% [H,P]=ttest2(MS_inRS.age_baseline,remaining_MS_fromRS.age_baseline)
% 
% nanmean(cell2mat(MS_inRS.disease_duration))
% nanmean(cell2mat(remaining_MS_fromRS.disease_duration))
% [H,P]=ttest2(cell2mat(MS_inRS.disease_duration),cell2mat(remaining_MS_fromRS.disease_duration))
% 
% mode((MS_inRS.gender))
% mode((remaining_MS_fromRS.gender))
% [H,P]=ttest2((MS_inRS.gender),(remaining_MS_fromRS.gender))
% 
% %% comparison HA in AD e viceversa
% 
% mean(RS_inAD.age)
% mean(remaining_RS_fromAD.age)
% [H,P]=ttest2(RS_inAD.age,remaining_RS_fromAD.age)
% 
% mode(RS_inAD.sex)
% mode(remaining_RS_fromAD.sex)
% [H,P]=ttest2(RS_inAD.sex,remaining_RS_fromAD.sex)
% 
% nanmean(RS_inAD.mmse_at_scan)
% nanmean(remaining_RS_fromAD.mmse_at_scan)
% [H,P]=ttest2(RS_inAD.mmse_at_scan,remaining_RS_fromAD.mmse_at_scan)
% 
% mode(RS_inAD.dementia_incident_type)
% mode(remaining_RS_fromAD.dementia_incident_type)
% [H,P]=ttest2(RS_inAD.dementia_incident_type,remaining_RS_fromAD.dementia_incident_type)
% 
% mode(RS_inAD.apoe4)
% mode(remaining_RS_fromAD.apoe4)
% [H,P]=ttest2(RS_inAD.apoe4,remaining_RS_fromAD.apoe4)
% 
% %% comparisons HA in MS vs HA not in MS
% 
% mean(RS_inMS.age)
% mean(remaining_RS_fromMS.age)
% [H,P]=ttest2(RS_inMS.age,remaining_RS_fromMS.age)
% 
% mode(RS_inMS.sex)
% mode(remaining_RS_fromMS.sex)
% [H,P]=ttest2(RS_inMS.sex,remaining_RS_fromMS.sex)
% 
% mean(RS_inMS.mmse_at_scan)
% mean(remaining_RS_fromMS.mmse_at_scan)
% [H,P]=ttest2(RS_inMS.mmse_at_scan,remaining_RS_fromMS.mmse_at_scan)
% 
% mode(RS_inMS.apoe4)
% mode(remaining_RS_fromMS.apoe4)
% [H,P]=ttest2(RS_inMS.apoe4,remaining_RS_fromMS.apoe4)
% 
% mode(RS_inMS.dementia_incident_type)
% mode(remaining_RS_fromMS.dementia_incident_type)
% [H,P]=ttest2(RS_inMS.dementia_incident_type,remaining_RS_fromMS.dementia_incident_type)
% 
