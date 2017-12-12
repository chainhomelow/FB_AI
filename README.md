# FB_AI
# Code samples for the Facebook AI residency
# 12/12/2017

## eyetracking_timecourse.m is an in-house script I built to do timecourse analyses on eye-movement data outputted from an SR ## Research Eyelink 1000 so we could see how long eyes focused on specific stimuli throughout an experimental trial, binning ## that trial in certan amounts of time. We did not have a way to do this analysis before this script. I wrote it for fun over ## the course of a weekend.


## pca_nuisance_regressors_final.m is another in-house script I wrote for our fMRI data, this time for a specific project. 
## Here I was running a Principle Components Analysis on a set of functional scans where the data were nuisance regressors 
## from all the scans. The purpose was to find common nuisance regressors via PCA and then remove them from all scans.

## ae01.acquisitiom_exinction_preprocessing.sh is a BASH script outlining one of our preprocessing pipelines for imaging data ## although each of our pipelines is tailored to a specific project
