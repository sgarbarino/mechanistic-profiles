
%% plot the results via tsne and ternary plot
% I merged together the data so I can plot AD, MS and RS together +
% separate the full data from the ES data from the HC data so they can be
% plotted 

%% tsne and ternary bootstrap on all subjects with DPM info
clear
load bootstrap_data.mat
figure(1)
tsne_bootstrap
figure(2)
triangle

%% tsne endstage data only
clear
load ES_data.mat 
figure(3)
tsne_endstage
figure(4)
triangle_endstage

%% tsne HC subjects only
clear
load HC_data.mat
figure(5)
tsne_HC
figure(6)
triangle_HC
