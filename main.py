import sys
# sys.path.append( '../GP_progression_model_V2/src' ) ## path where the functions of GPPM are
sys.path.append( '/Users/saragarbarino/Documents/WORK/MIDA/labs/liscomp/code/GPPM/GP_progression_model_V2/src' )
import GP_progression_model
import numpy as np
import torch
from torch.autograd import Variable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from copy import deepcopy
import bct
import functions_tp as tp

###### Load GPPM results (model_.pth; tr_.pth, misc.dat)
## Results are on a set of brain regions (NOTE: in this specific (toy) example are symmetrised, so BCT will be symmetrised too)

torch.backends.cudnn.benchmark = True
torch.backends.cudnn.fastest = True
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = GP_progression_model.GP_Progression_Model()
model.Load('./GPPM_data/')
# model.Plot()

## Load and assign GPs and time-ranges
tr = model.Return_time_parameters()
new_x = model.Transform_subjects()

x_min = float('inf')
x_max = float('-inf')
list_biomarker = model.names_biomarkers
for bio_pos, biomarker in enumerate(list_biomarker):
    bio_id = np.where([model.names_biomarkers[i] == biomarker for i in range(model.N_biomarkers)])[0][0]
    x_data = (model.model.time_reparameterization(model.x_torch)[bio_id]).detach().data.numpy()
    if np.float(np.min(x_data)) < x_min:
        x_min = np.float(np.min(x_data))
    if np.float(np.max(x_data)) > x_max:
        x_max = np.float(np.max(x_data))
x_range = Variable(torch.arange(x_min, x_max, np.float((x_max - x_min) / 150)))
x_range = x_range.reshape(x_range.size()[0], 1)

## predictions contain N_real realisations of each GP and their derivatives.
N_real = 200
predictions = [[] for i in range(len(list_biomarker))]
for bio_pos, biomarker in enumerate(list_biomarker):
    bio_id = np.where([model.names_biomarkers[i] == biomarker for i in range(model.N_biomarkers)])[0][0]
    predictions[bio_id] = []
    for i in range(N_real):
        predictions[bio_id].append(model.model.branches[bio_id][1](x_range))

###### Compute epicentres
# individual epicentres are the region with the highest value of disease severity in the last follow-up;
# the cohort-level epicenter is then the region that most frequently appears as individual-level epicenters.

Y_last = np.zeros((len(list_biomarker)-7, model.N_subs)) # NOTE: only on cortical regions, so if GPPM run on the whole brain, removing the subcortical regions
for biom in range(len(list_biomarker)-7):
    for sub in range(model.N_subs):
        Y_last[biom, sub] = model.y[biom][sub][np.argmax(model.x[biom][sub])]
counts = np.bincount(np.argmax(Y_last, 0))
epi = np.argmax(counts)
epi_name = model.names_biomarkers[epi]

####### Load connectome and cortical distances data (on cortical regions)
cortical = pd.DataFrame(pd.read_csv('./connectome_data/cortical_distances_original.csv', header=0, index_col=0)) # NOTE: The cortical distance cross-hemisphere is set to 0
connectome = pd.DataFrame(pd.read_csv('./connectome_data/connectome_mean.csv', header=0,index_col=0))
N = np.shape(connectome)[0]

## Normalisation of the connectome and the cortical distances to [0,1]
cortical_norm = cortical/np.max(cortical)
connectome_norm = connectome/np.max(connectome)

####### Graph-theory metrics computation, symmetrisation cross-hemisphere and normalisation in [0,1]
## Weighted degree
gt_deg = bct.strengths_und(connectome_norm.values)
gt_deg_sym = (gt_deg[:int(N/2)] + gt_deg[int(N/2):])/2

## Betweenness centrality
gt_betw = bct.betweenness_wei(1/connectome_norm.values)
gt_betw_sym = (gt_betw[:int(N/2)] + gt_betw[int(N/2):])/2

## Closeness centrality
gt_close = (len(connectome.values)-1)/np.sum(bct.distance_wei(1./connectome_norm.values)[0],1)
gt_close_sym = (gt_close[:int(N/2)] + gt_close[int(N/2):])/2

## Clustering coefficient
gt_clust = bct.clustering_coef_wu(connectome_norm.values)
gt_clust_sym = (gt_clust[:int(N/2)] + gt_clust[int(N/2):])/2

## Inverse of weighted degree
gt_invdeg = 1./gt_deg
gt_invdeg_sym = (gt_invdeg[:int(N/2)] + gt_invdeg[int(N/2):])/2

## Inverse of clustering coefficient
gt_invclust = 1./gt_clust
gt_invclust_sym = (gt_invclust[:int(N/2)] + gt_invclust[int(N/2):])/2

## Shortest path to epicenter (connectome-based)
# NOTE: connectome distance is computed after taking the mean over the two hemispheres
connectome_norm_sym = (connectome_norm.values[int(N/2):,int(N/2):] + connectome_norm.values[:int(N/2), :int(N/2)])/2
dist_conn_sym = bct.distance_wei(1./connectome_norm_sym);
gt_short_conn_sym = 1 - dist_conn_sym[0][epi,:] # the closest = the most vulnerable
# if not, uncomment the following and comment the previous:
# dist_conn = bct.distance_wei(1./connectome_norm.values)
# gt_short_conn=dist_conn[0][epi,:]
# gt_short_conn_sym = (gt_short_conn[:int(N/2)] + gt_short_conn[int(N/2):])/2

## Shortest path to epicenter (cortical-based)
# NOTE: cortical distance cannot be computed cross-hemisphere, so taking the mean over the two hemispheres
cortical_norm_sym = (cortical_norm.values[int(N/2):,int(N/2):] + cortical_norm.values[:int(N/2), :int(N/2)])/2
dist_cort_sym = bct.distance_wei(1./cortical_norm_sym);
gt_short_cort_sym = 1 - dist_cort_sym[0][epi,:] # the closest = the most vulnerable

## Constant progression
gt_constant_sym = np.ones((int(N/2),))

## Normalisation in [0,1]

gt_deg_sym = tp.norm(gt_deg_sym)
gt_betw_sym = tp.norm(gt_betw_sym)
gt_close_sym = tp.norm(gt_close_sym)
gt_clust_sym = tp.norm(gt_clust_sym)
gt_invdeg_sym = tp.norm(gt_invdeg_sym)
gt_invclust_sym = tp.norm(gt_invclust_sym)
gt_short_conn_sym = tp.norm(gt_short_conn_sym)
gt_short_cort_sym = tp.norm(gt_short_cort_sym)

####### Topological profiles
## Matrix sigma
sigma = np.vstack([gt_deg_sym, gt_betw_sym, gt_close_sym, gt_clust_sym, gt_invdeg_sym, gt_invclust_sym,
                   gt_short_conn_sym, gt_short_cort_sym, gt_constant_sym]).T
gt_metrics = ['weighted degree','betweenness','closeness','clustering','inverse degree', 'inverse clustering',
              'shortest path (connectome)','shortest path (cortical)','constant']

#### Cohort-level profiles
## Vector Y of maximum of the cohort-level trajectories (on cortical regions)
predictions = predictions[:int(N/2)]

## Parameters for ISRA
beta0 = np.ones((np.shape(sigma)[1],1))
maxiter = 10000

# Run the algorithm for cohort-level profiles (for each realisation)
beta = [[] for i in range(N_real)]
for i in range(N_real):
    dGP = [predictions[bio][i][:, 1].detach().numpy() for bio in range(len(list_biomarker)-7)]
    GP = [predictions[bio][i][:, 0].detach().numpy() for bio in range(len(list_biomarker)-7)]
    Y = np.expand_dims(np.max(dGP, 1), axis=1)
    beta[i] = tp.isra(sigma, Y, beta0, maxiter);
    print(i)

for i in range(200):
    beta[i][[beta[i]<1e-8]] = 0.
print("Cohort-level profiles\n")
for i in range(len(beta[0])):
    print("{first} = {second:.3f} +/- {third:2f}".format(first=gt_metrics[i], second=np.mean(beta,0)[i].squeeze(), third=np.std(beta,0)[i].squeeze()))

plt.figure()
plt.plot(range(len(Y)),Y)
plt.errorbar(range(len(Y)),np.matmul(sigma,np.mean(beta,0)).squeeze(),yerr=np.matmul(sigma,np.std(beta,0)).squeeze(), color='red')
plt.savefig('./output_cohort.png')

#### Individual profiles
## Vector of individual Y values
Y_indiv = np.zeros((model.N_biomarkers-7,model.N_subs))
for sub in range(model.N_subs):
    x_data = new_x[0][sub] # Individual temporal data (same for each marker if there's no missing data)
    # compute the time at which the (average) individual time intersects dGP
    time_inter_indiv = np.argmin(np.abs(x_range.numpy() * model.x_mean_std[0][1] + model.x_mean_std[0][0] - np.mean(
        x_data * model.x_mean_std[0][1] + model.x_mean_std[0][0])))
    for bio in range(model.N_biomarkers-7):
        dGP_mean = np.mean([predictions[bio][i][:, 1].detach().numpy() for i in range(N_real)], 0)
        Y_indiv[bio,sub] = dGP_mean[time_inter_indiv]

# Run the algorithm for individual profiles (using same parameters as cohort-level)
beta_indiv = [[] for i in range(model.N_subs)]
for sub in range(model.N_subs):
    beta_indiv[sub] = tp.isra(sigma, np.expand_dims(Y_indiv[:,sub], axis=1), beta0, maxiter)
    print(sub)

for i in range(model.N_subs):
    beta_indiv[i][[beta_indiv[i]<1e-8]] = 0.
print("\nAverage individual-level profiles\n")
for i in range(len(beta[0])):
    print("{first} = {second:.3f} +/- {third:2f}".format(first=gt_metrics[i], second=np.mean(beta_indiv,0)[i].squeeze(), third=np.std(beta_indiv,0)[i].squeeze()))


###### Visualisation using TSNE

from sklearn.manifold import TSNE

centers = np.mean(beta,0) # average cohort-level profile
X = np.array(beta_indiv).squeeze()
X = np.vstack((X,centers.T)) # stacked together

X_tnse = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(X)

plt.figure()
plt.scatter(X_tnse[:,0],X_tnse[:,1])
plt.scatter(X_tnse[-1,0],X_tnse[-1,1], color='red')
plt.savefig('./tnse.png')
