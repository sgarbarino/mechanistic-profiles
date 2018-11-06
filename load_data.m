%% loading step preliminar 

if strcmp(dataset,'ADNI')
    disp('Do not have access to raw AD data!')
    %     clean_ADNI; % load and pre-process ADNI data
    %     select_ROIs_ADNI; % AD: select all ROIs (i.e. h=1)
    %     zscoresROIs = (table2array(zscoresPAT_ROIs)'); % separate HC and PAT
    %     zscoresROIsHC = (table2array(zscoresHC_ROIs)');
    %     zscoresROIsALL = [zscoresROIs zscoresROIsHC];
    %     %% clean and load connectome
    %     script_mean_load_connectome_FS;
    %     script_mean_clean_connectome;
    %     clear idx p h
    %     %% merge all the information in 13 macro-regions:
    %     merge_lobes;
    %     zscoresROIs = zscores_volumes_NEW_lobes';
    %     connect = connect_lobes;
    %     ROIs = connectome_label_lobes;
    %     connectome_label = connectome_label_lobes';
    %     %     connectome_mean_lobes = mean(connect_lobes,4);
    %     connectome_std_lobes = std(connect_lobes,[],4);
    %     generate_GP_input;
    
elseif strcmp(dataset,'RS')
    disp('Do not have access to raw RS data!')
    %     clean_RS; % load and pre-process HA data
    %     Data = Data_pure; Data_demo = Data_demo_pure; Data_vol = Data_value_pure; % select just subjects with no cortical infarct + dementia at BL
    %     zscores_volumes_NEW_table = zscores_volumes_NEW_table_pure;zscores_volumes_NEW = zscores_volumes_NEW_pure;volumesall_NEW = volumes_NEW_pure;
    %     select_ROIs_RS; % HA:  select all ROIs (i.e. h=1)
    %     zscoresROIs = (table2array(zscores_ROIs)');
    %     %% clean and load connectome
    %     script_mean_load_connectome_FS;
    %     script_mean_clean_connectome;
    %     clear idx p h
    %     %% merge all the information in 13 macro-regions:
    %     merge_lobes;
    %     zscoresROIs = zscores_volumes_NEW_lobes';
    %     connect = connect_lobes;
    %     ROIs = connectome_label_lobes;
    %     connectome_label = connectome_label_lobes';
    %     %     connectome_mean_lobes = mean(connect_lobes,4);
    %     connectome_std_lobes = std(connect_lobes,[],4);
    %     generate_GP_input;
    
elseif strcmp(dataset,'PPMS')
    disp('Do not have access to raw PPMS data!')
    %     clean_PPMS; % load and pre-process PPMS data
    %     select_ROIs_PPMS; % MS:  select all ROIs (i.e. h=1)
    %     zscoresROIs = (table2array(zscoresPAT_ROIs)'); % for all
    %     zscoresROIsHC = (table2array(zscoresHC_ROIs)');
    %     %% clean and load connectome
    %     script_mean_load_connectome_FS;
    %     script_mean_clean_connectome;
    %     clear idx p h
    %     %% merge all the information in 13 macro-regions:
    %     merge_lobes;
    %     zscoresROIs = zscores_volumes_NEW_lobes';
    %     connect = connect_lobes;
    %     ROIs = connectome_label_lobes;
    %     connectome_label = connectome_label_lobes';
    %     %     connectome_mean_lobes = mean(connect_lobes,4);
    %     connectome_std_lobes = std(connect_lobes,[],4);
    %     generate_GP_input;
    
elseif strcmp(dataset,'dummy')
    clean_dummy;
    select_ROIs_dummy; % dummy: select all ROIs (i.e. h=1)
    zscoresROIs = (table2array(zscoresPAT_ROIs)'); % separate HC and PAT
    zscoresROIsHC = (table2array(zscoresHC_ROIs)');
    zscoresROIsALL = [zscoresROIs zscoresROIsHC];
    %% clean and load connectome
    script_mean_load_connectome_FS;
    script_mean_clean_connectome;
    clear idx p h
    %% merge all the information in 13 macro-regions:
    merge_lobes;
    zscoresROIs = zscores_volumes_NEW_lobes';
    connect = connect_lobes;
    ROIs = connectome_label_lobes;
    connectome_label = connectome_label_lobes';
    %     connectome_mean_lobes = mean(connect_lobes,4);
    connectome_std_lobes = std(connect_lobes,[],4);
    generate_GP_input;

else
    print('Error! Dataset not recognised')
    
end


