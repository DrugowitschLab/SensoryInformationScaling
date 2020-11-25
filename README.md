# Information Scaling in Large Neural Populations

This repository contains the scripts to reproduce the figures in the main text (Fig 2- 8) and supplementary information (Fig S1-S16) of the following article:

MohammadMehdi Kafashan, Anna Jaffe, Selmaan N Chettih, Ramon Nogueira, I&ntilde;igo Arandia-Romero, Christopher D Harvey, Rub&eacute;n Moreno-Bote, Jan Drugowitsch (2020). [Scaling of sensory information in large neural populations reveals signatures of information-limiting correlations]. Under review.


## Preliminaries

All scripts have been tested and run under MATLAB R2020a on MacOS Catalina. They should also work on other operating systems or later MATLAB versions. They might also run on [GNU Octave](https://www.gnu.org/software/octave/), the free MATLAB alternative, but this has not been tested.

The provided figure-generating scripts don't run _as is_. They require per-computing Fisher information moments, fit parameters, and other interim results that take up `>200Gb`, such that we don't provide them here. We provide all the scripts required to compute these interim results. See _Usage_ for how to compute them.

## Installation

All scripts are written in MATLAB. To run the scripts, first, download a [copy of the scripts](https://github.com/DrugowitschLab/SensoryInformationScaling/archive/main.zip), and extract them to a folder of your choice. Furthermore, download the [neural data](https://doi.org/10.6084/m9.figshare.13274951) from Figshare and place it into the `data` folder. 


## Usage

The main directory has several folders including `data`, `figs`, `fits`, `moment_cache`, `reordering_cache`, `shared`, `simData` and `tuning_fits` folders. Make sure that the MATLAB path points to the folder to which you extracted the scripts to. Add the `shared` folder to the MATLAB path as there are multiple utility functions in this folder that are used to perform the analysis. 

In the manuscript we number mice from 1 to 6. These numbers map to the mouse identifiers used in the scripts as follows:
```
mouse 1: m25  (2 sessions, a-b; 10% contrast)
mouse 2: m26  (2 sessions, a-b; 10% contrast)
mouse 3: aj42 (5 sessions, a-e; 10% contrast)
mouse 4: aj43 (7 sessions, a-g; 10% contrast)
mouse 5: aj60 (4 sessions, a-d; 10% & 25% contrast)
mouse 6: aj61 (3 sessions, a-c; 10% & 25% contrast)
```
Each identifier is followed by a letter that indicates the session. For example, the data for the second session of mouse 3 has the identifier `aj42b`.

## Computing interim results

Before executing the scripts that generate the figure plots, we need to generate a bunch of interim results that these scripts rely on. These results need to be computed for all mice/sessions and pairs of orientations. Computing these can take a considerable amount of time, and the use of computing clusters (if available) is recommended.

### Tuning curve fits

To perform tuning fits for recorded neural data we use the `singleNeurAnlys` function. For instance, to perform tuning curve fits for mouse/session `m25b`, we need to execute the following command in the MATLAB workspace:
```
singleNeurAnlys('m25b');
```

### Fisher information moments
First, we need to compute the first and second-order moments of Fisher information for different population sizes. We can do so for mouse/session `m25b` by
```
ori1id     = 1;         %   index of the first orientation in the virtual discrimination task
ori2id     = 2;         %   index of the second orientation in the virtual discrimination task
conid      = 1;         %   contrast index
popsamples = 10000;     %   number of bootstrap samples
computeMoments('m25b', ori1id, ori2id, conid, popsamples);
```
In the above,`ori1id`, `ori2id` index of the first/second orientation while `conid` represents the contrast index (for a dataset with _X_ different contrast `conid` ranges from 1 to _X_). The above script need to be re-run for all mice/sessions, all possible orientation pairs (by changing `ori1id` and `ori2id`), and contrasts (by changing `conid`). 

The `computeMoments` function bootstraps information increase estimates, and writes the bootstrap samples to the `moment_cache` folder (to a file named `m25b_o1-2_c1.mat` for the above example). These bootstrap samples form the basis for fitting the information increase curves.

To perform the analysis for trial-shuffled data, across trials, we need to provide a value for the `subsample` argument of the `computeMoments` script:
```
subsample   = 'shuf';
computeMoments('m25b', ori1id, ori2id, conid, popsamples, subsample);
```
Above command writes the bootstrap samples to a file named `m25b_o1-2_c1_shuf.mat` in the `moment_cache` folder. 

To perform the analysis only for trials with high running speed we need to instead call `computeMoments` as follows:
```
subsample   = 'hispd';
computeMoments('m25b', ori1id, ori2id, conid, popsamples, subsample);
```
The same analysis for low running speed-trials is performed by setting `subsample` to `lospd`.

### Information scaling model fits

Once the Fisher information moments have been computed and saved in `moment_cache`, we use the `scalemodel_fit` script to draw posterior samples for various information scaling models for a single dataset and single pair of information as follows:
```
scalemodel_fit('m25b_o1-2_c1.mat');
```
This command fits the information scaling for mouse/session `m25b`, orientation identifiers 1 and 2 (`o1_2`), and contrast identifier 1 (`c1`). Specific subsets of the data (e.g., trials with high running speed), or trial-shuffled data can be fit by calling the above script for the respective moment cache file (e.g., `m25b_o1-2_c1_shuf.mat`).

To fit information scaling models to multiple orientation pairs simultanously, we use the script `scalemodel_fit_multi`. The first argument to this function specifies the name of the joint analysis (which will be the used filename to store the fits), and the remaining arguments specify the moment files that the fits are performed for. For example,
```
scalemodel_fit_multi("m25b_dori1_c1", 'm25b_o1-2_c1','m25b_o1-8_c1','m25b_o2-3_c1','m25b_o3-4_c1','m25b_o4-5_c1','m25b_o5-6_c1','m25b_o6-7_c1','m25b_o7-8_c1')
```
performs the fit across 8 different orientation pairs, all separated by `45°`, for mouse/session `m25b`.

### Information estimates for optimal reordering of neurons
To compute greedy optimal orderings of neurons in terms of the information we utilize `computeOrderings` function for a given dataset, discrimination tasks, and contrast. The function accept a boolean argument named `shuf` to perform the analysis for trial-shuffled data. For instance, to perform the information estimates for optimal reordering of neurons for dataset `m25` and discrimination task between direction `45°` and `90°` we need to execute the following commands in the MATLAB workspace:
```
ori1id      = 1;        %   index of the first orientation in the virtual discrimination task
ori2id      = 2;        %   index of the second orientation in the virtual discrimination task
conid       = 1;        %   contrast index
computeOrderings('m25b', ori1id, ori2id, conid);
```

To perform the same computation for trial-shuffled data we need to execute the following commands:
```
ori1id      = 1;        %   index of the first orientation in the virtual discrimination task
ori2id      = 2;        %   index of the second orientation in the virtual discrimination task
conid       = 1;        %   contrast index
shuffle     = true;     %   perform trial-shuffling of neural data before optimal re-ordering analysis
computeOrderings('m25b', ori1id, ori2id, conid, shuffle);
```

Results of the above analysis will be stored in `reordering_cache`. This analysis needs to be repeated for different mice/sessions, orientation pairs, and contrasts. 

## Main text and supplementary figures

### Figures in the main text
For reproducing the figures in the main text execute the following commands in the workspace:
```
plotFig2;
plotFig3;
plotFig4;
plotFig6;
plotFig7;
plotFig8;
```

### Figures in the supplementary information (SI)
Except for the figures that are explicitly discussed below, all other SI figures should be executed by typing `plotFigS#`, where `#` needs to be replaced by the SI figure number. 


### Fig S7
We first need to generate the simulated neural activity using a Gaussian population model, and Linear-nonlinear-Poisson (LNP) population model ([Kanitscheider et al. (2015)](https://doi.org/10.1073/pnas.1508738112)). We utilize `simPopActivity` function to generate the simulated neural activity and the neural activity of each simulation is written to `simData` folder with the name of `sim[simid]` where `simid` is the integer representing the simulation id and ranges from 1 to 4. Simulation details for different `simid` numbers are expressed as follows:
- `simid = 1` - Gaussian population with limited information, $$N=1000$$, $$T=1000$$, $$I_\infty=20$$
- `simid = 2` - Gaussian population with unlimited information,  $$N=1000$$, $$T=1000$$, $$I_\infty=\infty$$
- `simid = 3` - LNP population with limited information, $$N=2500$$, $$T=2500$$
- `simid = 4` - LNP population with unlimited information, $$N=2500$$, $$T=2500$$

See Sec. 4.1 and Sec 4.2 of the SI for more details on the multivariate Gaussian population model and LNP model. After storing the simulated neural activity, `computeMoments` is utilized to compute the Fisher information moments for the simulated data. For instance, to compute the information moments for `sim1`, using only `N=300` neurons and `T=500` trials, we utilized the following command:

```
N           = 300;      %   number of neurons to consider
T           = 500;      %   number of trials to consider
popsamples  = 10000;    %   number of bootstrap samples
computeMoments('sim1', N, T, popsamples);
```
After that, we perform the information scaling estimates using `scalemodel_fit` function and then execute the `plotFigS7` in the MATLAB workspace to generate this SI figure.

### Figure S12
First need to execute the single neuron analysis for all sessions of mice 5 and 6 using the `singleNeurAnlys('aj60a')` function as follows:
```
subsample = 'none';
singleNeurAnlys_mouse5n6('aj60a', subsample);
```
After that, we execute `plotFigS12` in the MATLAB workspace to generate this SI figure.

### Figure S15
First execute the function `sim_eye_mvm` to simulate neural population activity, model of [Kanitscheider et al. (2015)](https://doi.org/10.1073/pnas.1508738112), with `N=1000` and `T=20000` and different combination of `pn` and `pt` as follows:
```
sim_eye_mvm
```
This function will store the information increase moments to the `moment_cache`. We then utilize the `scalemodel_fit` function to the stored moment data in the `moment_cache` folder to estimate the information scaling for this simulated data. After that, we execute `plotFigS15` in the MATLAB workspace to generate this SI figure.

### Figure S16
First execute the function `simIcov` to simulate populations of different sizes `M` as described in Sec. 4.1 of the SI. For that we need to execute the following command in MATLAB workspace for `M=500, 1000, 5000, 10000, 15000`:
```
M   = 500;  %   Population size of the simulated neural activity. 
simIcov(M)
```
The covariance statistics will be written to `cov_cache` folder. After executing the above command for differnet `M`, we execute `plotFigS16` in the MATLAB workspace to generate this SI figure.