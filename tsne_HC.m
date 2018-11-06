%% tsne with HC subjects only

AD_pp_HC = [nodAD troAD traAD proAD]';
MS_pp_HC = [nodMS troMS traMS proMS]';
RS_pp = [nodRS troRS traRS proRS]'; % for RS  it's the entirety of subjects

% %%
meas_HC = [AD_pp_HC'; MS_pp_HC'; RS_pp']; 
for i = 1:length(AD_pp_HC)
    labels{i} = 'AD-HC';
end
for i = length(AD_pp_HC)+1:length(AD_pp_HC)+length(MS_pp_HC)
    labels{i} = 'PPMS-HC';
end
for i = length(AD_pp_HC)+length(MS_pp_HC)+1:length(AD_pp_HC)+length(MS_pp_HC)+length(RS_pp)
    labels{i} = 'HA';
end

labels = labels';

meas_HC_cent = [meas_HC ; centerAD; [0,0,0,0]; centerRS];
labels_cen_HC{1} = 'AD cohort profile';
labels_cen_HC{2} = 'PPMS cohort profile';
labels_cen_HC{3} = 'HA cohort profile';
labels_HC_cent = [labels; labels_cen_HC'];

tic;
rng('default') % for reproducibility
Y_euc_all_HC = tsne(meas_HC_cent,'Algorithm','barneshut','Distance','euclidean','Perplexity',140,'Verbose',2);
toc;

gscatter(Y_euc_all_HC(:,1),Y_euc_all_HC(:,2),labels_HC_cent,['r','g','b','r','g','b'],...
    ['.','.','.','p','p','p'],[5,5,5,25,25,25]);

