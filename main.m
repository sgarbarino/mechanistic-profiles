%% script to run the whole pipeline

% STEPS 1-4 can be run with dummy dataset only, which mimics ADNI data.

% STEPS 5-9 can be run on dummy, or AD, HA or PPMS processed,
% unrecognisable data. Such data can be requested to the author at:
% sara.garbarino@inria.fr

% STEP 10 run with AD+HA+PPMS data to generate the plots; it does not run
% with dummy dataset, as it is a single dataset and ternary and tsne would
% not be meaningful.

%% 1: load and pre-process AD/HA/PPMS data; 
% load and pre-process connectomes and compute GT metrics

% They run on raw data, if there's access to them. 
% To test them, I created a dummy dataset that mimics AD data

clear all

% dataset = 'dummy';
dataset = 'ADNI';
% dataset = 'PPMS';
% dataset = 'RS';

load_data; % load data + connectomes

if strcmp(dataset, 'dummy')
    ranking; % create GT metrics
end

%% 2: GP progression modelling step - http://gpprogressionmodel.inria.fr/ 
% requires a CSV to be uploaded with biomarkers (for us z-scores of 13 volumes), 
% diagnosis, and age

%% 3: Load the .npy data and put them in matlab-friendly format. 

% the output of the GP progression are group-level trajectories + individual 
% information on the time shift. 
% Finally data are saved in a .mat, cleaned of raw data information, for
% sharing (see below)
if strcmp(dataset, 'dummy')
    load_GP
end
%% 4: run the "endstage" part
if strcmp(dataset, 'dummy') 
    load_ES
end

%% 5: load ADNI, RS or PPMS data without raw data but with GP output - for sharing
if strcmp(dataset, 'ADNI')
    load data/AD.mat 
elseif strcmp(dataset, 'PPMS')
    load data/PPMS.mat
elseif strcmp(dataset, 'RS')
    load data/HA.mat 
end

%% 6: run mechanistic profiles model Y_meas_cont = beta*X_all_cont

% X contains all the GT metrics: 
% Betwenness centrality (SB); closeness centrality (SCLO); weighted degree
% (SD); clustering coefficient (SCLUT); inverse of degree (SID); inverse of
% clustering coefficient (SCLU); shortest path (SSH); proximity spread
% (SPR); constant propagation (ones).
GT_metrics ={'betwenness', 'closeness', 'degree','clustering','inverse degree','inverse clust',...
    'shortest','proximity','constant'};

X_all_cont = [SB.mean SCLO.mean SD.mean SCLUT.mean...
    SID.mean SCLU.mean ...
    SSH.mean SPR.mean ...
    ones(size(SPR.mean))];

% Y_meas_cont and _meas_cont contain the values at population-level (mean 
% and std), already normalized

Y_meas_cont = table_OUT.der_values;
Y_meas_cont_std = table_OUT.der_values_std;

% ISRA
beta0_a = ones(size(X_all_cont,2),size(Y_meas_cont,2));
opt.maxiter = 50000;

mean_beta_a_cont = isra(X_all_cont, X_all_cont', Y_meas_cont, beta0_a, opt);
std_beta_a_cont =  isra(X_all_cont, X_all_cont', Y_meas_cont_std, beta0_a, opt);

formatSpec = 'Weights for GT metric %s: mean (std) %8.4f  (%8.4f)\n';

for i = 1:9
    fprintf(formatSpec,GT_metrics{i}, mean_beta_a_cont(i),std_beta_a_cont(i));
end

% residuals
Y_output_isra = (mean_beta_a_cont'*X_all_cont')';
residuals_isra = sum((Y_meas_cont - Y_output_isra).^2);

% generate N_boot Y_meas_cont data from data distribution

N_boot = 100;
opt.maxiter = 10000;
% [bootstat,bootsam]=bootstrp(N_boot,@mean,Y_meas_cont); % bootsam encodes the index of the subjects

for i = 1:N_boot
    % select database
%     Y_meas_cont_boot_cont = Y_meas_cont(bootsam(:,i),:);   
    Y_meas_cont_boot_cont = Y_meas_cont_std .*randn(1) + Y_meas_cont;
    % re-run ISRA
    beta_a_cont_boot_cont(:,i) = isra(X_all_cont, X_all_cont', Y_meas_cont_boot_cont, beta0_a, opt);
    disp(i)
end
mean_beta_a_cont_boot_cont = mean(beta_a_cont_boot_cont,2);
std_beta_a_cont_boot_cont = std(beta_a_cont_boot_cont')';

formatSpec = 'Weights for GT metric %s: mean (std) %8.4f  (%8.4f)\n';
for i = 1:9
    fprintf(formatSpec,GT_metrics{i}, mean_beta_a_cont_boot_cont(i),std_beta_a_cont_boot_cont(i));
end

% permutation testing for each beta
for i = 1:9
    [aa(i),bb(i),cc(i)] = permutationTest(beta_a_cont_boot_cont(i,:),zeros(1,100),...
    10000,'plotresult', 1, 'showprogress', 250);
end

%% 7: run individual mechanistic profiles and individual bootstrap
for ind = 1:length(runs)
    Y_meas_cont_id(:,ind) = struct_OUT_all.(runs{ind}).der_values ;
end

beta0_a_ind = 5*ones(size(X_all_cont,2),size(Y_meas_cont_id,2));
opt.maxiter = 10000;

beta_a_cont_ind = isra(X_all_cont, X_all_cont', Y_meas_cont_id, beta0_a_ind, opt);

N_boot = 100;

[bootstat,bootsam]=bootstrp(N_boot,@mean,Y_meas_cont_id); % bootsam encodes the index of the subjects

beta_a_cont_boot = zeros([size(beta_a_cont_ind),N_boot]);
mean_beta_a_cont_boot = zeros(size(beta_a_cont_ind,1),N_boot);
opt.maxiter = 10000;

for i = 1:N_boot
    % select database
    Y_meas_cont_boot = Y_meas_cont_id(bootsam(:,i),:);
    % re-run ISRA
    beta_a_cont_boot(:,:,i) = isra(X_all_cont, X_all_cont', Y_meas_cont_boot, beta0_a_ind, opt);
    mean_beta_a_cont_boot(:,i) = mean(beta_a_cont_boot(:,:,i),2); 
    disp(i)
end
mean_mean_beta_a_cont_boot = mean(beta_a_cont_ind,2);
std_beta_a_cont_boot = std(beta_a_cont_ind')';

formatSpec = 'Weights for GT metric %s: mean (std) %8.4f  (%8.4f)\n';
for i = 1:9
    fprintf(formatSpec,GT_metrics{i}, mean_mean_beta_a_cont_boot(i),std_beta_a_cont_boot(i));
end

%% 8: run mechanistic profiles model for end-stage data

Y_meas_cont_ES = table_OUT_es.YY; % mean of all subj in struct_OUT_ES
Y_meas_cont_std = table_OUT_es.Y_STD;

opt.maxiter = 10000;

% population-level
mean_beta_a_es = isra(X_all_cont, X_all_cont', Y_meas_cont_ES, beta0_a, opt);
std_beta_a_es =  isra(X_all_cont, X_all_cont', Y_meas_cont_std, beta0_a, opt);

formatSpec = 'Weights for GT metric %s: mean (std) %8.4f  (%8.4f)\n';
for i = 1:9
    fprintf(formatSpec,GT_metrics{i}, mean_beta_a_es(i),std_beta_a_es(i));
end

% subject-wise
for ind = 1:length(runs)
    Y_meas_cont_id_es(:,ind) = struct_OUT_ES{ind}.subj ;
end

beta0_a_ind = 5*ones(size(X_all_cont,2),size(Y_meas_cont_id_es,2));
opt.maxiter = 100000;

beta_a_cont_ind_es = isra(X_all_cont, X_all_cont', Y_meas_cont_id_es, beta0_a_ind, opt);


%% 9: randomization of graphs 
%1. compute GT metrics at random (i.e. on random N-nodes graph)
load data/connectomes.mat
G=WattsStrogatz(N,floor(N/2)-2,1);
adj_rand = full(adjacency(G));
ranking_rnd;

%2. compute mechanisms correlating atrophy on subjects and random GT
X_all_cont_rnd = [SB_rnd.mean SCLO_rnd.mean SD_rnd.mean SCLUT_rnd.mean...
    SID_rnd.mean SCLU.mean ...
    SSH_rnd.mean SPR_rnd.mean ...
    ones(size(SPR_rnd.mean))];
mean_beta_a_cont_rnd = isra(X_all_cont_rnd, X_all_cont_rnd', Y_meas_cont, beta0_a, opt);
std_beta_a_cont_rnd =  isra(X_all_cont_rnd, X_all_cont_rnd', Y_meas_cont_std, beta0_a, opt);

formatSpec = 'Weights for randomized GT metric %s: mean (std) %8.4f  (%8.4f)\n';
for i = 1:9
    fprintf(formatSpec,GT_metrics{i}, mean_beta_a_cont_rnd(i),std_beta_a_cont_rnd(i));
end
Y_output_isra_rnd = (mean_beta_a_cont'*X_all_cont_rnd')';
residuals_isra_rnd = sum((Y_meas_cont - Y_output_isra_rnd).^2);

%3. obtain individual mechanisms
for ind = 1:length(runs)
    Y_meas_cont_id(:,ind) = struct_OUT_all.(runs{ind}).der_values ;
end

beta_a_cont_ind_rnd = isra(X_all_cont_rnd, X_all_cont_rnd', Y_meas_cont_id, beta0_a_ind, opt);

beta_a_cont_boot_rnd = zeros([size(beta_a_cont_ind_rnd),N_boot]);
mean_beta_a_cont_boot_rnd = zeros(size(beta_a_cont_ind_rnd,1),N_boot);
opt.maxiter = 10000;

for i = 1:N_boot
    Y_meas_cont_boot = Y_meas_cont_id(bootsam(:,i),:);
    beta_a_cont_boot_rnd(:,:,i) = isra(X_all_cont_rnd, X_all_cont_rnd', Y_meas_cont_boot, beta0_a_ind, opt);
    mean_beta_a_cont_boot_rnd(:,i) = mean(beta_a_cont_boot_rnd(:,:,i),2); 
    disp(i)
end
mean_mean_beta_a_cont_boot_rnd = mean(beta_a_cont_ind_rnd,2);
std_beta_a_cont_boot_rnd = std(beta_a_cont_ind_rnd')';

formatSpec = 'Weights for GT metric %s: mean (std) %8.4f  (%8.4f)\n';
for i = 1:9
    fprintf(formatSpec,GT_metrics{i}, mean_mean_beta_a_cont_boot_rnd(i),std_beta_a_cont_boot_rnd(i));
end


%% compute centers and population for plots - saved in separate mat files for endstage, HC or all subjects - only for dummy
% centerAD = [sum(mean_beta_a_cont_rnd(1:4)); sum(mean_beta_a_cont_rnd(5:6)); mean_beta_a_cont_rnd(7); mean_beta_a_cont_rnd(8)];
% centerMS = [sum(mean_beta_a_cont_rnd(1:4)); sum(mean_beta_a_cont_rnd(5:6)); mean_beta_a_cont_rnd(7); mean_beta_a_cont_rnd(8)];
% centerRS = [sum(mean_beta_a_cont_rnd(1:4)); sum(mean_beta_a_cont_rnd(5:6)); mean_beta_a_cont_rnd(7); mean_beta_a_cont_rnd(8)];
% nodAD = sum(beta_a_cont_ind_rnd(1:4,:));
% troAD = sum(beta_a_cont_ind_rnd(5:6,:));
% traAD = beta_a_cont_ind_rnd(7,:);
% proAD = beta_a_cont_ind_rnd(8,:);
% nodMS = sum(beta_a_cont_ind_rnd(1:4,:));
% troMS = sum(beta_a_cont_ind_rnd(5:6,:));
% traMS = beta_a_cont_ind_rnd(7,:);
% proMS = beta_a_cont_ind_rnd(8,:);
% nodRS = sum(beta_a_cont_ind_rnd(1:4,:));
% troRS = sum(beta_a_cont_ind_rnd(5:6,:));
% traRS = beta_a_cont_ind_rnd(7,:);
% proRS = beta_a_cont_ind_rnd(8,:);
% centerAD_boot = [sum(beta_a_cont_boot_rnd(1:4)); ...
%     sum(beta_a_cont_boot_rnd(5:6)); beta_a_cont_boot_rnd(7);...
%     beta_a_cont_boot_rnd(8)];
% centerPPMS_boot = [sum(beta_a_cont_boot_rnd(1:4)); ...
%     sum(beta_a_cont_boot_rnd(5:6)); beta_a_cont_boot_rnd(7);...
%     beta_a_cont_boot_rnd(8)];
% centerRS_boot = [sum(beta_a_cont_boot_rnd(1:4)); ...
%     sum(beta_a_cont_boot_rnd(5:6)); beta_a_cont_boot_rnd(7);...
%     beta_a_cont_boot_rnd(8)];
% meanAD_boot = [sum(mean_mean_beta_a_cont_boot_rnd(1:4)); ...
%     sum(mean_mean_beta_a_cont_boot_rnd(5:6)); mean_mean_beta_a_cont_boot_rnd(7);...
%     mean_mean_beta_a_cont_boot_rnd(8)];
% stdAD_boot = [sum(beta_a_cont_ind_rnd(1:4)); ...
%     sum(beta_a_cont_ind_rnd(5:6)); beta_a_cont_ind_rnd(7);...
%     beta_a_cont_ind_rnd(8)];
% meanPPMS_boot = [sum(mean_mean_beta_a_cont_boot_rnd(1:4)); ...
%     sum(mean_mean_beta_a_cont_boot_rnd(5:6)); mean_mean_beta_a_cont_boot_rnd(7);...
%     mean_mean_beta_a_cont_boot_rnd(8)];
% stdPPMS_boot = [sum(beta_a_cont_ind_rnd(1:4)); ...
%     sum(beta_a_cont_ind_rnd(5:6)); beta_a_cont_ind_rnd(7);...
%     beta_a_cont_ind_rnd(8)];
% meanRS_boot = [sum(mean_mean_beta_a_cont_boot_rnd(1:4)); ...
%     sum(mean_mean_beta_a_cont_boot_rnd(5:6)); mean_mean_beta_a_cont_boot_rnd(7);...
%     mean_mean_beta_a_cont_boot_rnd(8)];
% stdRS_boot = [sum(beta_a_cont_ind_rnd(1:4)); ...
%     sum(beta_a_cont_ind_rnd(5:6)); beta_a_cont_ind_rnd(7);...
%     beta_a_cont_ind_rnd(8)];
% centerAD_ES = [sum(mean_beta_a_es(1:4)); sum(mean_beta_a_es(5:6)); mean_beta_a_es(7); mean_beta_a_es(8)];
% centerMS_ES = [sum(mean_beta_a_es(1:4)); sum(mean_beta_a_es(5:6)); mean_beta_a_es(7); mean_beta_a_es(8)];
% centerRS_ES = [sum(mean_beta_a_es(1:4)); sum(mean_beta_a_es(5:6)); mean_beta_a_es(7); mean_beta_a_es(8)];
% meanAD_es = [sum(beta_a_cont_ind_es(1:4)); ...
%     sum(beta_a_cont_ind_es(5:6)); beta_a_cont_ind_es(7);...
%     beta_a_cont_ind_es(8)];
% meanPPMS_es = [sum(beta_a_cont_ind_es(1:4)); ...
%     sum(beta_a_cont_ind_es(5:6)); beta_a_cont_ind_es(7);...
%     beta_a_cont_ind_es(8)];
% meanRS_es = [sum(beta_a_cont_ind_es(1:4)); ...
%     sum(beta_a_cont_ind_es(5:6)); beta_a_cont_ind_es(7);...
%     beta_a_cont_ind_es(8)];

%% 10: plot for tsne and ternary

plot_data;


