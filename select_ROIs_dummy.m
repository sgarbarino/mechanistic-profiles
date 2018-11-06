%% select ROIs (all of them)
 
demoAD = demoPAT(demoPAT.DX =='Dementia'|demoPAT.DX =='MCI to Dementia'|demoPAT.DX =='NL to Dementia',:);
h = ones(size(volumesall_NEW,2),1);
% ROIs selection:
ROIs = volumesall_NEW.Properties.VariableNames(find(h)); % names of the selected ROIs


%% select the corresponding zscores for both PAT and HC

for j = 1: length(ROIs)
    zscoresPAT_ROIs(:,j) = zscores_volumesPAT_NEW_table(:,strcmp(zscores_volumesPAT_NEW_table.Properties.VariableNames, ROIs{j})); %here I can decide if using PAT or one subtype (they must be computed, as SPMS in line 11)
    zscoresHC_ROIs(:,j) = zscores_volumesHC_NEW_table(:,strcmp(zscores_volumesHC_NEW_table.Properties.VariableNames, ROIs{j})); 
end
%%
zscoresPAT_ROIs.Properties.VariableNames = ROIs;
zscoresHC_ROIs.Properties.VariableNames = ROIs;


clear j 

clear  zscores_volumesHC_NEW_table zscores_volumesLMCI_NEW_table zscores_volumesAD_NEW_table zscores_volumesPAT_NEW_table