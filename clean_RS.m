%% load FS Volumes on RS subjects + demo info
%% clean the data: select just Volumes, merge L/R + over lobes

% load RS dataset

load data_RS_FS.mat
load names_volumes_FS.mat 

%% here select if thickness or volumes + WM parts

idx_vol = strfind(Data.Properties.VariableNames,'volume');
for i = 1:length(idx_vol)
    if isempty(idx_vol{1,i})
        idx_vol{1,i}=0;
    end
end
clear i

idx_WM = zeros(size(idx_vol));
idx_WM([9:12,16,17,19,27:33])=1;

idx_std = strfind(Data.Properties.VariableNames,'thicknessstd');  % have to remove thicknessstd
for i = 1:length(idx_std)
    if isempty(idx_std{1,i})
        idx_std{1,i}=0;
    end
end

idx_thick = strfind(Data.Properties.VariableNames,'thickness'); 
for i = 1:length(idx_thick)
    if isempty(idx_thick{1,i})
        idx_thick{1,i}=0;
    end
end
clear i

idx_thick = cell2mat(idx_thick)+(idx_WM)-cell2mat(idx_std);
idx_vol = cell2mat(idx_vol)+(idx_WM);

Data_tot_vol = Data(:,logical((idx_vol)));
Data_tot_thick = Data(:,logical(idx_thick));
clear i idx_thick idx_thick1 idx_thick2 idx_vol

%% select if vol or thick
Data_tot_value = Data_tot_vol; 

clear Data_tot_vol Data_tot_thick
Data_tot_demo = Data(:,[1:3,59,402:408]);

% change names
Data_tot_value.Properties.VariableNames = Names_Volumes;
clear Names_Volumes
% remove what is not present in connectome_label and reorder ROIs to match
% connectome label
order = [15:48,1:14,49:82];
Data_tot_value = Data_tot_value(:,order);
clear order

%% remove cortical infarct and dementia_prevalent
Data_tot_value(Data_tot_demo.corticalinfarct>0,:)=[];
Data_tot_demo(Data_tot_demo.corticalinfarct>0,:)=[];
Data_tot_value(Data_tot_demo.dementia_prevalent>0,:)=[];
Data_tot_demo(Data_tot_demo.dementia_prevalent>0,:)=[];

%% re-combine Data merging vol and demo
clear Data
Data_tot = [Data_tot_demo Data_tot_value];

%% regress out also TIV, gender

X = [ones(size(Data_tot.IntraCranialVol)) Data_tot.IntraCranialVol  Data_tot.sex];

beta = pinv(X)*Data_tot_value{:,:};
Y_teo = X*beta;
volumes_NEW_tot = array2table(Data_tot_value{:,:} - Y_teo);% + mean_YHC_teo);
volumes_NEW_tot.Properties.VariableNames = Data_tot_value.Properties.VariableNames;

%% now compute z scores of grey matter volumes

mean_all = mean(volumes_NEW_tot{:,:});
std_all = std(volumes_NEW_tot{:,:});

% zscores
zscores_volumes_NEW_tot = (volumes_NEW_tot{:,:} - ...
    repmat(mean_all, size(volumes_NEW_tot{:,:},1),1))./...
    repmat(std_all, size(volumes_NEW_tot{:,:},1),1);
zscores_volumes_NEW_tot_table = array2table(zscores_volumes_NEW_tot);
zscores_volumes_NEW_tot_table.Properties.VariableNames = Data_tot_value.Properties.VariableNames;

clear  Y_teo mean_Y_teo X beta std_all mean_all
%% now repeat removing dementia_incident (just pure aging)
    
Data_pure = Data_tot;
Data_pure(Data_tot_demo.dementia_incident2015>0,:)=[];
Data_demo_pure = Data_pure(:,1:11);
Data_value_pure = Data_pure(:,12:end);

X = [ones(size(Data_pure.IntraCranialVol)) Data_pure.IntraCranialVol  Data_pure.sex];
beta = pinv(X)*Data_value_pure{:,:};
Y_teo = X*beta;
volumes_NEW_pure = array2table(Data_value_pure{:,:} - Y_teo);% + mean_YHC_teo);
volumes_NEW_pure.Properties.VariableNames = Data_value_pure.Properties.VariableNames;
mean_all = mean(volumes_NEW_pure{:,:});
std_all = std(volumes_NEW_pure{:,:});
zscores_volumes_NEW_pure = (volumes_NEW_pure{:,:} - ...
    repmat(mean_all, size(volumes_NEW_pure{:,:},1),1))./...
    repmat(std_all, size(volumes_NEW_pure{:,:},1),1);
zscores_volumes_NEW_table_pure = array2table(zscores_volumes_NEW_pure);
zscores_volumes_NEW_table_pure.Properties.VariableNames = Data_value_pure.Properties.VariableNames;


clear  Y_teo mean_Y_teo X beta std_all mean_all




