%% select ROIs (all of them)

demoOLD = Data_demo(Data_demo.age > 74,:);
demoYOUNG = Data_demo(Data_demo.age < 54,:);

zscores_volumesOLD_NEW_table = zscores_volumes_NEW_table(Data_demo.age > 74,:);
zscores_volumesYOUNG_NEW_table = zscores_volumes_NEW_table(Data_demo.age < 54,:);

h = ones(size(volumesall_NEW,2),1);
ROIs = Data_vol.Properties.VariableNames(find(h)); % names of the selected ROIs

%% select the corresponding zscores for both PAT and HC

for j = 1: length(ROIs)
    zscores_ROIs(:,j) = zscores_volumes_NEW_table(:,strcmp(zscores_volumes_NEW_table.Properties.VariableNames, ROIs{j})); 
end
%%
zscores_ROIs.Properties.VariableNames = ROIs;


clear j 
clear zscores_volumesOLD_NEW_table zscores_volumesYOUNG_NEW_table demoYOUNG zscores_volumes_NEW