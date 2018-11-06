%% script to compute GT metrics

addpath('/home/sgarbari/Documents/WORK/POND/code/BCT/2016_01_16_BCT/')
%% graph-related measure

% nodal stress = highest degree
degree = sum(connectome,2);

% tropic failure = lowest degree
inversedegree = 1./degree;

% clustering
C_true = clustering_coef_wu(connectome);


% INVERSE clustering  
C = 1./clustering_coef_wu(connectome);

% betweenness centrality 
%(Node betweenness centrality is the fraction of all shortest paths in the 
%network that contain a given node. Nodes with high values of betweenness 
%centrality participate in a large number of shortest paths.)
B = betweenness_wei(1./connectome); 

% closeness centrality 
dist = distance_wei(1./connectome);

L = (length(connectome)-1)./(sum(dist,2));
connectome_ordinal = [1:length(connectome_label)]';
N = length(connectome_label);%height(nodal_stress_sorted);
char = sum(dist,2);


%% tables 

nodal_stress = table(degree,connectome_label,connectome_ordinal);
tropic_failure = table(inversedegree,connectome_label, connectome_ordinal);
clust = table(C,connectome_label,connectome_ordinal);
clust_true = table(C_true,connectome_label,connectome_ordinal);
betw = table(B,connectome_label,connectome_ordinal);
clos = table(L,connectome_label,connectome_ordinal);
char_length = table(char,connectome_label,connectome_ordinal);


nodal_stress_sorted = nodal_stress;
tropic_failure_sorted = tropic_failure;
clust_sorted = clust;
clust_true_sorted = clust_true;
betw_sorted = betw;
clos_sorted = clos;
char_length_sorted  = char_length;

