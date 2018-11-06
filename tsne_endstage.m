%% tsne end-stage

AD_pp_ES = [nodAD_ES troAD_ES traAD_ES proAD_ES]';
MS_pp_ES = [nodMS_ES troMS_ES traMS_ES proMS_ES]';
RS_pp_ES = [nodRS_ES troRS_ES traRS_ES proRS_ES]';

% %%
meas_ES = [AD_pp_ES'; MS_pp_ES'; RS_pp_ES']; 

meas_ES_all = [meas_ES ; centerAD_ES; centerMS_ES; centerRS_ES];
labels_cen_ES{1} = 'AD end-stage cohort profile';
labels_cen_ES{2} = 'PPMS end-stage cohort profile';
labels_cen_ES{3} = 'HA end-stage cohort profile';
labels_ES_all = [labels; labels_cen_ES'];

for i = 1:length(AD_pp_ES)
    labels{i} = 'AD-ES';
end
for i = length(AD_pp_ES)+1:length(AD_pp_ES)+length(MS_pp_ES)
    labels{i} = 'PPMS-ES';
end
for i = length(AD_pp_ES)+length(MS_pp_ES)+1:length(AD_pp_ES)+length(MS_pp_ES)+length(RS_pp_ES)
    labels{i} = 'HA-ES';
end

labels = labels';
opt.MaxIter=5000;
tic;
rng('default') % for reproducibility
Y_euc_all_ES = tsne(meas_ES_all,'Algorithm','barneshut','Distance','euclidean','Perplexity',100,'Verbose',2);
toc;

gscatter(Y_euc_all_ES(:,1),Y_euc_all_ES(:,2),labels_ES_all,['r','g','b','r','g','b'],...
    ['.','.','.','p','p','p'],[5,5,5,25,25,25]);