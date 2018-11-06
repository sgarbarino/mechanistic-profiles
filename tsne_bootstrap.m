%% tsne all subject + bootstrap

AD_pp = [nodAD troAD traAD proAD]';
MS_pp = [nodMS troMS traMS proMS]';
RS_pp = [nodRS troRS traRS proRS]';

meas = [AD_pp'; MS_pp'; RS_pp']; 

for i = 1:length(AD_pp)
    labels{i} = 'AD individual profiles';
end
for i = length(AD_pp)+1:length(AD_pp)+length(MS_pp)
    labels{i} = 'PPMS individual profiles';
end
for i = length(AD_pp)+length(MS_pp)+1:length(AD_pp)+length(MS_pp)+length(RS_pp)
    labels{i} = 'HA individual profiles';
end

meas_boot = [centerAD_boot  ; centerPPMS_boot  ; centerHA_boot ; centerAD ; centerMS ; centerRS]; 

for i = 1:length(centerAD_boot)
    labels_boot{i} = 'AD bootstrapped cohort profiles';
end
for i = length(centerAD_boot)+1:length(centerAD_boot)+length(centerPPMS_boot)
    labels_boot{i} = 'PPMS bootstrapped cohort profiles';
end
for i = length(centerAD_boot)+length(centerPPMS_boot)+1:length(centerAD_boot)+length(centerPPMS_boot)+length(centerHA_boot)
    labels_boot{i} = 'HA bootstrapped cohort profiles';
end

labels_boot{i+1} = 'AD cohort profile';
labels_boot{i+2} = 'PPMS cohort profile';
labels_boot{i+3} = 'HA cohort profile';

meas_all = [meas ; meas_boot];
labels_all = [labels'; labels_boot'];
 
clear Y_euc_all
tic;
rng('default') % for reproducibility
Y_euc_all = tsne(meas_all,'Algorithm','barneshut','Distance','euclidean','Perplexity',100,'Verbose',2);
toc;


gscatter(Y_euc_all(:,1),Y_euc_all(:,2),labels_all,['r','g','b','r','g','b','r','g','b'],...
    ['.','.','.','p','p','p','p','p','p'],[5,5,5,15,15,15,25,25,25])