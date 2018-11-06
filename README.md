# mechanistic-profiles

Code for building mechanistic profiles of neurological diseases from grey matter volumetric information.

The code is in Matlab and tested with V 2017 a). It requires the "Brain Connectivity Toolbox" (https://sites.google.com/site/bctnet/), "npy-matlab" (https://github.com/kwikteam/npy-matlab), the "ternary plot" (https://fr.mathworks.com/matlabcentral/fileexchange/7210-ternary-plots) and tSNE (https://fr.mathworks.com/help/stats/tsne.html) Matlab adds-on to run.


## Usage

Code can be run from the main.m script; a dummy data set reproducing ADNI data "dummy_data.mat" is available to test the functioning of the code. For the GP Progression Model part, please refer to http://gpprogressionmodel.inria.fr/; you would need to upload a csv file containing information on the GM volumes, the diagnosis and the age, which is generated authomatically by the main.m.

