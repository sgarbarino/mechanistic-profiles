# topological-profiles

Code for building topological profiles of neurological diseases from brain regional information (from MRI or PET imaging)

***05/2023: major update of the code to run with the released version of GPPM***

## Note
* The topological profile code is now implemented in python (3.7) for a smoother interface with GPPM:  https://gitlab.inria.fr/epione/GP_progression_model_V2.
* It has a number of standard requirementes (numpy, pandas, torch, matplotlib) and bctpy: https://pypi.org/project/bctpy/.
* Also, it requires GPPM to be locally in your computer, as it uses some of its data loading procedures for loading GPPM data.

## Usage
* Code can be run from the main.py script; please note that it will require (standard) saved output from GPPM in the form of a "model_.pth" and "tr_.pth" files.
* By default it runs with some (toy) example of ADNI data, for which GPPM has produced a set of output which can be found in the "GPPM_data" folder
* By default it computes a set of 9 graph theory metrics on the structural connectome/cortical distance;
  * weighted degree
  * betweeness centrality
  * closeness centrality
  * clustering coefficient
  * inverse of weighted degree
  * inverse of clustering coefficient
  * shortest path to epicenter (via connections)
  * shortest path to epicenter (via cortical paths - also called "cortical proximity")
  * constant, homogeneous, isotropic propagation
* The folder "connectome_data" contains some pre-computed average structural connectome and average cortical distance data on a set of HCP brains: https://www.humanconnectomeproject.org/ (for how these were computed, see Oxtoby et al, Frontiers, 2017 and Garbarino et al, eLife 2019)

-------

## Older version (Matlab)
The previous version of the code, running on Matlab, can be still found in the "matlab" branch. Please note that it will run smoothly on *dummy data*, but IT WILL NOT EASILY RUN if you have data obtained with the released version GPPM: https://gitlab.inria.fr/epione/GP_progression_model_V2. Main.m launches the code on dummy ADNI data.

It required:
* "Brain Connectivity Toolbox" (https://sites.google.com/site/bctnet/);
* "npy-matlab" (https://github.com/kwikteam/npy-matlab)
