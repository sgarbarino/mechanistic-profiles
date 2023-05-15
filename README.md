# topological-profiles

Code for building mechanistic profiles of neurological diseases from grey matter volumetric information. 

***05/2023: major update of the code to run with the released version of GPPM: https://gitlab.inria.fr/epione/GP_progression_model_V2***

## Note
* The topological profile code is now implemented in python (3.7) for a smoother interface with GPPM.
* It has a number of standard requirementes (numpy, pandas, torch, matplotlib) and bctpy: https://pypi.org/project/bctpy/.

## Usage
Code can be run from the main.py script; please note that it will require (standard) saved output from GPPM in the form of a "model_.pth" and "tr_.pth" files.


## Older version (Matlab)
The previous version of the code, running on Matlab, can be still found in the "matlab" branch. Please note that it will run smoothly on *dummy data*, but IT WILL NOT EASILY RUN if you have data obtained with the released version GPPM: https://gitlab.inria.fr/epione/GP_progression_model_V2. Main.m launches the code on dummy ADNI data.

It required:
* "Brain Connectivity Toolbox" (https://sites.google.com/site/bctnet/);
* "npy-matlab" (https://github.com/kwikteam/npy-matlab)
