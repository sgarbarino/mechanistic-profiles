%% run ranking of connectomics measures 
addpath('BCT/')

%% inizialization

SALVA_BETW = zeros(size(connect,4),size(ROIs,2));
SALVA_CLOS = zeros(size(connect,4),size(ROIs,2));
SALVA_CLUST = zeros(size(connect,4),size(ROIs,2));
SALVA_CLUSTTRUE = zeros(size(connect,4),size(ROIs,2));
SALVA_DEGREE = zeros(size(connect,4),size(ROIs,2));
SALVA_INVDEGREE = zeros(size(connect,4),size(ROIs,2));
SALVA_SHORT = zeros(size(connect,4),size(ROIs,2));

SALVA_BETW_intensity = zeros(size(connect,4),size(ROIs,2));
SALVA_CLOS_intensity = zeros(size(connect,4),size(ROIs,2));
SALVA_CLUST_intensity = zeros(size(connect,4),size(ROIs,2));
SALVA_CLUSTTRUE_intensity = zeros(size(connect,4),size(ROIs,2));
SALVA_DEGREE_intensity = zeros(size(connect,4),size(ROIs,2));
SALVA_INVDEGREE_intensity = zeros(size(connect,4),size(ROIs,2));
SALVA_SHORT_intensity = zeros(size(connect,4),size(ROIs,2));

%% loop
tic
for j = 1: size(connect,4) % loop on graphs
    connectome = connect(:,:,:,j); 
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
    
    SALVA_BETW_intensity(j,:) = betw_sorted.B;
    SALVA_CLOS_intensity(j,:) = clos_sorted.L;
    SALVA_DEGREE_intensity(j,:) = nodal_stress_sorted.degree;
    SALVA_INVDEGREE_intensity(j,:) = tropic_failure_sorted.inversedegree;
    SALVA_CLUST_intensity(j,:) = clust_sorted.C;
    SALVA_CLUSTTRUE_intensity(j,:) = clust_true_sorted.C_true;
    SALVA_SHORT_intensity(j,:) = short_path_sorted.distances_epi;
    
end

disp('Completed on all the connectomes '); toc

%% now mean over subjects

SALVA_BETW_mean = mean(SALVA_BETW_intensity);SALVA_BETW_std = std(SALVA_BETW_intensity);
SB = table(SALVA_BETW_mean',SALVA_BETW_std',connectome_label);SB.Properties.VariableNames={'mean','std','connectome_label'};

SALVA_CLOS_mean = mean(SALVA_CLOS_intensity);SALVA_CLOS_std = std(SALVA_CLOS_intensity);
SCLO = table(SALVA_CLOS_mean',SALVA_CLOS_std',connectome_label);SCLO.Properties.VariableNames={'mean','std','connectome_label'};

SALVA_DEGREE_mean = mean(SALVA_DEGREE_intensity);SALVA_DEGREE_std = std(SALVA_DEGREE_intensity);
SD = table(SALVA_DEGREE_mean',SALVA_DEGREE_std',connectome_label);SD.Properties.VariableNames={'mean','std','connectome_label'};

SALVA_INVDEGREE_mean = mean(SALVA_INVDEGREE_intensity);SALVA_INVDEGREE_std = std(SALVA_INVDEGREE_intensity);
SID = table(SALVA_INVDEGREE_mean',SALVA_INVDEGREE_std',connectome_label);SID.Properties.VariableNames={'mean','std','connectome_label'};

SALVA_CLUST_mean = mean(SALVA_CLUST_intensity);SALVA_CLUST_std = std(SALVA_CLUST_intensity);
SCLU = table(SALVA_CLUST_mean',SALVA_CLUST_std',connectome_label);SCLU.Properties.VariableNames={'mean','std','connectome_label'};

SALVA_CLUSTTRUE_mean = mean(SALVA_CLUSTTRUE_intensity);SALVA_CLUSTTRUE_std = std(SALVA_CLUSTTRUE_intensity);
SCLUT = table(SALVA_CLUSTTRUE_mean',SALVA_CLUSTTRUE_std',connectome_label);SCLUT.Properties.VariableNames={'mean','std','connectome_label'};

SALVA_SHORT_mean = mean(SALVA_SHORT_intensity);SALVA_SHORT_std = std(SALVA_SHORT_intensity);
SSH = table(SALVA_SHORT_mean',SALVA_SHORT_std',connectome_label);SSH.Properties.VariableNames={'mean','std','connectome_label'};

%% now create proximity
dijkstra;
tmpiso = 0:1/12:1; % equally spaced probability
% now add std variation
[va,vb]=mode(SALVA_VOL);
stdiso = (1-1./(24./vb));
tmptab=table(new_vol_accazzo',tmpiso',stdiso',ROIs(new_vol_accazzo)');
tmptab=sortrows(tmptab);
SPR = tmptab(:,2:end);
SPR.Properties.VariableNames={'mean','std','connectome_label'};

%% now normalise all between 0 and 1

SB.mean = (SB.mean - min(SB.mean))./(max(SB.mean)-min(SB.mean));
SCLO.mean = (SCLO.mean - min(SCLO.mean))./(max(SCLO.mean)-min(SCLO.mean));
SD.mean = (SD.mean - min(SD.mean))./(max(SD.mean)-min(SD.mean));
SCLUT.mean = (SCLUT.mean - min(SCLUT.mean))./(max(SCLUT.mean)-min(SCLUT.mean));

SID.mean = (SID.mean - min(SID.mean))./(max(SID.mean)-min(SID.mean));
SCLU.mean = (SCLU.mean - min(SCLU.mean))./(max(SCLU.mean)-min(SCLU.mean));

SSH.mean = (SSH.mean - max(SSH.mean))./(min(SSH.mean)-max(SSH.mean));
SPR.mean = (SPR.mean - max(SPR.mean))./(min(SPR.mean)-max(SPR.mean));

