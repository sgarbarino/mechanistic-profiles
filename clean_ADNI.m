%% load FS Volumes on ADNI subjects and the spreadsheet with demo information.
%% clean the data: select just Volumes, merge L/R + over lobes

% load ADNI dataset

load data_ADNI_FS.mat
load names_volumes_FS.mat 

%% remove various bits 

% step 1 - clean to keep just info that will also be in the FS connectomes (84 ROIs).
% separate thicnkess and volume

idx_vol = strfind(Data.Properties.VariableNames,'Volume');
for i = 1:length(idx_vol)
    if isempty(idx_vol{1,i})
        idx_vol{1,i}=0;
    end
end 
clear i
Data_tot_vol = Data(:,logical(cell2mat(idx_vol)));
idx_thick1 = strfind(Data.Properties.VariableNames,'ThicknessAverage'); %thickness+WM pparcellation (volumetric)
for i = 1:length(idx_thick1)
    if isempty(idx_thick1{1,i})
        idx_thick1{1,i}=0;
    end
end
clear i
idx_thick2 = strfind(Data.Properties.VariableNames,'VolumeWM'); %thickness+WM pparcellation (volumetric)
for i = 1:length(idx_thick2)
    if isempty(idx_thick2{1,i})
        idx_thick2{1,i}=0;
    end
end
idx_thick = cell2mat(idx_thick1)+cell2mat(idx_thick2);
Data_tot_thick = Data(:,logical(idx_thick));
clear i idx_thick idx_thick1 idx_thick2 idx_vol
% remove subjects with 0s volumes (don't have FS) 

% here remove regions of cerebellum or outside of 84 regions in Names_Volumes
for i = 1:length(Names_Volumes)
    A_vol(i,:)=strfind(Data_tot_vol.Properties.VariableNames, Names_Volumes{i,1});
end
clear i
for i = 1:length(Names_Volumes)
    for j = 1:width(Data_tot_vol)
        if isempty(A_vol{i,j})
            A_vol{i,j}=0;
        end
    end
end
B_vol = cell2mat(A_vol);

for i = 1:length(Names_Volumes)
    A_thick(i,:)=strfind(Data_tot_thick.Properties.VariableNames, Names_Volumes{i,1});
end
clear i
for i = 1:length(Names_Volumes)
    for j = 1:width(Data_tot_thick)
        if isempty(A_thick{i,j})
            A_thick{i,j}=0;
        end
    end
end
B_thick = cell2mat(A_thick);
clear i j A_thick A_vol

vol_kept = sum(B_vol(:,:)>1); 
thick_kept = sum(B_thick(:,:)>1); 

Data_tot_vol = Data_tot_vol(:,logical(vol_kept));

Data_tot_value = Data_tot_vol; %%%% choose here if volumes or thickness

clear Data_tot_vol Data_tot_thick
Data_tot_demo = Data(:,[1:30,37:50,57:64,433:438]); % separate demo

% change names
Data_tot_value.Properties.VariableNames = Names_Volumes;
clear Names_volumes
% eorder ROIs to match order in connectome label
Data_tot_value(:,[37,78])=[];
order = [26	27	28	30	31	33	35	36	37	38	39	40	41	42	45	44	46	47	48	49	50	51	52	53	55	56	57	58	59	60	32	61	63	23	62	29	54	43	34	24	20	21	69	12	1	74	65	64	66	67	68	70	71	73	75	76	77	78	79	80	81	82	3	2	4	5	6	7	8	9	10	11	13	14	15	16	17	18	72	19	22	25];
Data_tot_value = Data_tot_value(:,order);
clear order

%% re-combine Data merging vol and demo and removing subj which have NaN ICV, or ICV too small (i.e. <10)
clear Data_tot
Data_tot = [Data_tot_demo Data_tot_value];
clear Data_tot_demo Data_tot_value
% replace AGE with true AGE - AGE+Years_bl
Data_tot.AGE = Data_tot.AGE+Data_tot.Years_bl;

% clear from subj with NaN ICV + the minimum
Data_tot(isnan(Data_tot.ICV),:)=[];
Data_tot(Data_tot.ICV==min(Data_tot.ICV),:)=[];

% replace male with 1 and female with 0
for i = 1 : height(Data_tot)
    if strcmp(Data_tot.PTGENDER(i),'Female')
        Data_tot.PTGENDER{i} = 0;
    end
    if strcmp(Data_tot.PTGENDER(i),'Male')
        Data_tot.PTGENDER{i} = 1;
    end
end
clear i
Data_tot.PTGENDER = cell2mat(Data_tot.PTGENDER);

% NaN viscode in 0
Data_tot.VISCODE(isnan(Data_tot.VISCODE)) = 0;
% disease type into categorical array
Data_tot.DX_bl = categorical(Data_tot.DX_bl);
Data_tot.DX = categorical(Data_tot.DX);

% separate again
demo = Data_tot(:,1:58);% ICV in demo
volumes = Data_tot(:,59:end);
clear ans 

%% regress out age, PTGENDER and ICV - use HC to find the linear fit parameter

demoHC = demo(demo.DX =='NL'| demo.DX =='MCI to NL',:);
volumesHC = volumes(demo.DX =='NL'| demo.DX =='MCI to NL',:);

X = [ones(size(demoHC.ICV)) demoHC.ICV demoHC.AGE demoHC.PTGENDER];

beta = pinv(X)*volumesHC{:,:};
YHC_teo = X*beta;
mean_YHC_teo = repmat(mean(YHC_teo),size(YHC_teo,1),1);

volumesHC_NEW = array2table(volumesHC{:,:} - YHC_teo);% + mean_YHC_teo);
volumesHC_NEW.Properties.VariableNames = volumes.Properties.VariableNames;

%% use these parameters to build new Patients data
demoPAT = demo(demo.DX =='Dementia' | demo.DX =='Dementia to MCI' | demo.DX =='MCI' | demo.DX =='MCI to Dementia'|demo.DX =='NL to Dementia' | demo.DX =='NL to MCI',:);
volumesPAT = volumes(demo.DX =='Dementia' | demo.DX =='Dementia to MCI' | demo.DX =='MCI' | demo.DX =='MCI to Dementia'|demo.DX =='NL to Dementia' | demo.DX =='NL to MCI',:);

X_PAT = [ones(size(demoPAT.ICV)) demoPAT.ICV demoPAT.AGE demoPAT.PTGENDER];
Y_teo = X_PAT*beta;
mean_Y_teo = repmat(mean(Y_teo),size(Y_teo,1),1);

volumesPAT_NEW = array2table(volumesPAT{:,:} - Y_teo);% + mean_Y_teo);
volumesPAT_NEW.Properties.VariableNames = volumes.Properties.VariableNames;

%% now compute z scores of grey matter volumes against ALL, using new data

volumesall_NEW = array2table([volumesHC_NEW{:,:} ;volumesPAT_NEW{:,:}]); % merge volumes HC and PAT

volumesall_NEW.Properties.VariableNames = volumes.Properties.VariableNames;
mean_all = nanmean(volumesall_NEW{:,:});
std_all = nanstd(volumesall_NEW{:,:});

% zscores of PAT
zscores_volumesPAT_NEW = (volumesPAT_NEW{:,:} - ...
    repmat(mean_all, size(volumesPAT_NEW{:,:},1),1))./...
    repmat(std_all, size(volumesPAT_NEW{:,:},1),1);
zscores_volumesPAT_NEW_table = array2table(zscores_volumesPAT_NEW);
zscores_volumesPAT_NEW_table.Properties.VariableNames = volumes.Properties.VariableNames;

% compute also zscores for HC
zscores_volumesHC_NEW = (volumesHC_NEW{:,:} - ...
    repmat(mean_all, size(volumesHC_NEW{:,:},1),1))./...
    repmat(std_all, size(volumesHC_NEW{:,:},1),1);
zscores_volumesHC_NEW_table = array2table(zscores_volumesHC_NEW);
zscores_volumesHC_NEW_table.Properties.VariableNames = volumes.Properties.VariableNames;

%%
clear X_PAT Y_teo mean_Y_teo mean_YHC_teo X beta YHC_teo std_all mean_all
clear volumes volumesHC volumesPAT zscores_volumesHC_NEW zscores_volumesPAT_NEW
