%% run ranking of connectomics measures 
addpath('BCT/')

%% inizialization

%% loop

connectome = adj_rand;
connectome = weight_conversion(connectome,'normalize'); % connectome normalization
compute_GT; % compute graph-related metrics

distances = distance_wei(1./connectome); % for shortest path (oss: with epi = most atrophied region)

for i = 1:size(zscoresROIs,2)
    [~,epi(i)]=min(zscoresROIs(:,i));
end
[epi_mean] = mode(epi);
distances_epi = distances(epi_mean,:)';
short_path = table(distances_epi,connectome_label,connectome_ordinal);
short_path_sorted = short_path;

SALVA_BETW_intensity_rnd = betw_sorted.B;
SALVA_CLOS_intensity_rnd = clos_sorted.L;
SALVA_DEGREE_intensity_rnd = nodal_stress_sorted.degree;
SALVA_INVDEGREE_intensity_rnd = tropic_failure_sorted.inversedegree;
SALVA_CLUST_intensity_rnd = clust_sorted.C;
SALVA_CLUSTTRUE_intensity_rnd = clust_true_sorted.C_true;
SALVA_SHORT_intensity_rnd = short_path_sorted.distances_epi;

%% now mean over subjects

SB_rnd = table(SALVA_BETW_intensity_rnd,connectome_label);SB_rnd.Properties.VariableNames={'mean','connectome_label'};

SCLO_rnd = table(SALVA_CLOS_intensity_rnd,connectome_label);SCLO_rnd.Properties.VariableNames={'mean','connectome_label'};

SD_rnd = table(SALVA_DEGREE_intensity_rnd,connectome_label);SD_rnd.Properties.VariableNames={'mean','connectome_label'};

SID_rnd = table(SALVA_INVDEGREE_intensity_rnd,connectome_label);SID_rnd.Properties.VariableNames={'mean','connectome_label'};

SCLU_rnd = table(SALVA_CLUST_intensity_rnd,connectome_label);SCLU_rnd.Properties.VariableNames={'mean','connectome_label'};

SCLUT_rnd = table(SALVA_CLUSTTRUE_intensity_rnd,SALVA_CLUSTTRUE_std',connectome_label);SCLUT_rnd.Properties.VariableNames={'mean','std','connectome_label'};

SSH_rnd = table(SALVA_SHORT_intensity_rnd,connectome_label);SSH_rnd.Properties.VariableNames={'mean','connectome_label'};

%% now create proximity
dijkstra;
tmpiso = 0:1/12:1; % equally spaced probability
% now add std variation
[va,vb]=mode(SALVA_VOL);
stdiso = (1-1./(24./vb));
tmptab=table(new_vol_accazzo',tmpiso',stdiso',ROIs(new_vol_accazzo)');
tmptab=sortrows(tmptab);
SPR_rnd = tmptab(:,[2,4]);
SPR_rnd.Properties.VariableNames={'mean','connectome_label'};

%% now normalise all between 0 and 1

SB_rnd.mean = (SB_rnd.mean - min(SB_rnd.mean))./(max(SB_rnd.mean)-min(SB_rnd.mean));
SCLO_rnd.mean = (SCLO_rnd.mean - min(SCLO_rnd.mean))./(max(SCLO_rnd.mean)-min(SCLO_rnd.mean));
SD_rnd.mean = (SD_rnd.mean - min(SD_rnd.mean))./(max(SD_rnd.mean)-min(SD_rnd.mean));
SCLUT_rnd.mean = (SCLUT_rnd.mean - min(SCLUT_rnd.mean))./(max(SCLUT_rnd.mean)-min(SCLUT_rnd.mean));

SID_rnd.mean = (SID_rnd.mean - min(SID_rnd.mean))./(max(SID_rnd.mean)-min(SID_rnd.mean));
SCLU_rnd.mean = (SCLU_rnd.mean - min(SCLU_rnd.mean))./(max(SCLU_rnd.mean)-min(SCLU_rnd.mean));

SSH_rnd.mean = (SSH_rnd.mean - max(SSH_rnd.mean))./(min(SSH_rnd.mean)-max(SSH_rnd.mean));
SPR_rnd.mean = (SPR_rnd.mean - max(SPR_rnd.mean))./(min(SPR_rnd.mean)-max(SPR_rnd.mean));

