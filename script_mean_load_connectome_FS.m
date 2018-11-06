%% load connectome with idx away:

idx =[35,84]; %cerebellum 
%% load connectome
[connectome_mean,connectome_std,connectome_label,connect] = load_connectomeHCP_meanFS(idx);
