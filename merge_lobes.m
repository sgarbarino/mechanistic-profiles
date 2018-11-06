
% I want:

% Frontal 1
% Temporal 2
% Parietal 3
% Occipital 4
% Cingulate 5
% Insula 6
% Thalamus 7
% Caudate 8
% Putamen 9
% Pallidum 10
% Hippocampus 11
% Amygdala 12
% Accumbens 13

% I have to map ROIs to those areas - have to clean connectome AND
% zscoresROIs

size_lobes = 13;
    
    point = [2,5,1,4,2,2,3,2,5,4,1,4,1,2,2,1,1,1,1,4,3,5,1,3,5,1,1,3,2,3,1,2,2,6,7,8,9,10,11,12,13]';
    connect_lobes = zeros(size_lobes,size_lobes,1,size(connect,4));
    
    for conn_i = 1:size(connect,4)
        adj = zeros(size_lobes,size_lobes);
        for ii = 1:size_lobes
            for jj = (ii+1):size_lobes
                adj(ii,jj) = sum(sum(connect(point==ii,point==jj,1,conn_i)));
            end
        end
        adj = adj + adj.' - 2*diag(diag(adj));
        connect_lobes(:,:,1,conn_i) = adj;
    end
    clear ii jj conn_i
    
    connectome_label_lobes = {'Frontal','Temporal','Parietal','Occipital','Cingulate'...
        'Insula','Thalamus','Caudate','Putamen','Pallidum','Hippocampus','Amygdala','Accumbens'};
    connectome_mean_lobes = mean(connect_lobes,4);
    connectome_std_lobes = std(connect_lobes,[],4);
    
    if exist('volumesHC_NEW') % in AD      
        zscores_volumesHC_NEW_lobes = zeros(size(zscoresHC_ROIs,1),size_lobes);
        zscores_volumesPAT_NEW_lobes = zeros(size(zscoresPAT_ROIs,1),size_lobes);
        for i_vol = 1:size_lobes
            zscores_volumesHC_NEW_lobes(:,i_vol) = sum(zscoresHC_ROIs{:,point==i_vol},2)/size(zscoresHC_ROIs{:,point==i_vol},2);
            zscores_volumesPAT_NEW_lobes(:,i_vol) = sum(zscoresPAT_ROIs{:,point==i_vol},2)/size(zscoresPAT_ROIs{:,point==i_vol},2);
        end
        zscores_volumes_NEW_lobes = [zscores_volumesHC_NEW_lobes; zscores_volumesPAT_NEW_lobes];
        zscores_volumes_NEW_table_lobes = array2table(zscores_volumes_NEW_lobes);
        zscores_volumes_NEW_table_lobes.Properties.VariableNames = connectome_label_lobes;
        
    else % in HA
        zscores_volumes_NEW_lobes = zeros(size(zscores_volumes_NEW_pure,1),size_lobes);
        for i_vol = 1:size_lobes
            zscores_volumes_NEW_lobes(:,i_vol) = sum(zscores_volumes_NEW_table_pure{:,point==i_vol},2)/size(zscores_volumes_NEW_table_pure{:,point==i_vol},2);
        end
        zscores_volumes_NEW_table_lobes = array2table(zscores_volumes_NEW_lobes);
        zscores_volumes_NEW_table_lobes.Properties.VariableNames = connectome_label_lobes;
    end
    clear  i_vol std_all mean_all


