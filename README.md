# mechanistic-profiles

Code for building mechanistic profiles of neurological diseases from grey matter volumetric information. 

It runs on dummy data set mimiking ADNI data.

The code is in Matlab and tested with V 2017 a). 

## Requirements
* "Brain Connectivity Toolbox" (https://sites.google.com/site/bctnet/);
* "npy-matlab" (https://github.com/kwikteam/npy-matlab)

## Usage

Code can be run from the main.m script; a dummy data set mimiking ADNI data "dummy_data.mat" is available to test the functioning of the code. For the GP Progression Model part, please refer to http://gpprogressionmodel.inria.fr/; you would need to upload a csv file containing information on the GM volumes, the diagnosis and the age, which is generated authomatically by the main.m.

