%% load FS Volumes on PPMS subjects + demo info
%% clean the data: select just Volumes, merge L/R + over lobes

% load PPMS dataset

load data_PPMS_FS.mat
load names_volumes_FS.mat 

%% add AGE with true AGE - AGE+Years_bl
Data_tot.AGE = Data_tot.age_baseline+Data_tot.years;

% clear from subj with NaN ICV + the minimum
Data_tot(isnan(Data_tot.EstimatedTotalIntraCranialVol),:)=[];

% replace male with 1 and female with 0
for i = 1 : height(Data_tot)
    if strcmp(Data_tot.gender(i),'female')
        Data_tot.gender{i} = 0;
    end
    if strcmp(Data_tot.gender(i),'male')
        Data_tot.gender{i} = 1;
    end
end
clear i
Data_tot.gender = cell2mat(Data_tot.gender);

% NaN EDSS in 0
% Data_tot.EDSS(isnan(Data_tot.EDSS)) = 0;
% disease type into categorical array
Data_tot.subtype = categorical(Data_tot.subtype);


% separate again
demo = Data_tot(:,[1:10 ,end-1:end]);% ICV in demo
volumes = Data_tot(:,11:end-2);
clear ans 

%% regress out age, PTGENDER and ICV - use HC to find the linear fit parameter

demoHC = demo(demo.subtype =='HC',:);
volumesHC = volumes(demo.subtype =='HC',:);

X = [ones(size(demoHC.EstimatedTotalIntraCranialVol)) demoHC.EstimatedTotalIntraCranialVol demoHC.AGE demoHC.gender];

beta = pinv(X)*volumesHC{:,:};
YHC_teo = X*beta;
mean_YHC_teo = repmat(mean(YHC_teo),size(YHC_teo,1),1);

volumesHC_NEW = array2table(volumesHC{:,:} - YHC_teo);% + mean_YHC_teo);
volumesHC_NEW.Properties.VariableNames = volumes.Properties.VariableNames;

% use these parameters to build new Patients data
demoPAT = demo(demo.subtype =='PPMS' ,:);
volumesPAT = volumes(demo.subtype =='PPMS' ,:);

X_PAT = [ones(size(demoPAT.EstimatedTotalIntraCranialVol)) demoPAT.EstimatedTotalIntraCranialVol demoPAT.AGE demoPAT.gender];
Y_teo = X_PAT*beta;
mean_Y_teo = repmat(mean(Y_teo),size(Y_teo,1),1);

volumesPAT_NEW = array2table(volumesPAT{:,:} - Y_teo);% + mean_Y_teo);
volumesPAT_NEW.Properties.VariableNames = volumes.Properties.VariableNames;


%% now compute z scores of grey matter volumes against ALL (as Arman), using new data

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

clear X_PAT Y_teo mean_Y_teo mean_YHC_teo X beta YHC_teo std_all mean_all
clear volumes volumesHC volumesPAT zscores_volumesHC_NEW zscores_volumesPAT_NEW
