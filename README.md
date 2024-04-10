# iPSC-CM mathematical model calibration
This repository contains code to run simulations, model calibrations, and reproduce the figures in "Creating cell-specific computational models of stem cell-derived cardiomyocytes using optical experiments" by Yang J et al. This work was performed in the Cardiac Systems Pharmacology Lab of Dr. Eric A. Sobie in the Department of Pharmacological Sciences at the Icahn School of Medicine Mount Sinai.
## Abstract (in progress)
Human induced pluripotent stem cell-derived cardiomyocytes (iPSC-CMs) have gained traction as a powerful model in cardiac disease and therapeutics research, since iPSCs are self-renewing and can be derived from healthy and diseased patients without invasive surgery. However, current iPSC-CM differentiation methods produce cardiomyocytes with immature, fetal-like electrophysiological phenotypes, and the variety of maturation protocols in the literature results in phenotypic differences between labs. Heterogeneity of iPSC donor genetic backgrounds contributes to additional phenotypic variability. Several mathematical models of iPSC-CM electrophysiology have been developed to help understand the ionic underpinnings of, and to simulate, various cell responses, but these models individually do not capture the phenotypic variability observed in iPSC-CMs. Here, we tackle these limitations by developing a computational pipeline to calibrate cell preparation-specific iPSC-CM electrophysiological parameters. 

We used the genetic algorithm (GA), a heuristic parameter calibration method, to tune ion channel parameters in a mathematical model of iPSC-CM physiology. To systematically optimize an experimental protocol that generates sufficient data for parameter calibration, we created simulated datasets by applying various protocols to a population of in silico cells with known conductance variations, and we fitted to those datasets. We found that calibrating models to voltage and calcium transient data under 3 varied experimental conditions, including electrical pacing combined with ion channel blockade and changing buffer ion concentrations, improved model parameter estimates and model predictions of unseen channel block responses. This observation held regardless of whether the fitted data were normalized, suggesting that normalized fluorescence recordings, which are more accessible and higher throughput than patch clamp recordings, could sufficiently inform conductance parameters. Therefore, this computational pipeline can be applied to different iPSC-CM preparations to determine cell line-specific ion channel properties and understand the mechanisms behind variability in perturbation responses. 
## Code requirements
MATLAB version 2019 or higher; version 2019a was used to create the Kernik model populations, and 2020b-2023b were used to perform analysis.
## Models included within this repository:
Kernik DC, Morotti S, Wu H, Garg P, Duff HJ, Kurokawa J, et al. A computational model of induced pluripotent stem‐cell derived cardiomyocytes incorporating experimental variability from multiple data sources. J Physiol. 2019;597: 4533–4564. doi:10.1113/JP277724
O’Hara T, Virág L, Varró A, Rudy Y. Simulation of the Undiseased Human Cardiac Ventricular Action Potential: Model Formulation and Experimental Validation. PLoS Comput Biol. 2011;7: e1002061. doi:10.1371/journal.pcbi.1002061
Paci M, Hyttinen J, Aalto-Setälä K, Severi S. Computational models of ventricular- and atrial-like human induced pluripotent stem cell derived cardiomyocytes. Ann Biomed Eng. 2013;41: 2334–2348. doi:10.1007/s10439-013-0833-3
Paci M, Pölönen R-P, Cori D, Penttinen K, Aalto-Setälä K, Severi S, et al. Automatic Optimization of an in Silico Model of Human iPSC Derived Cardiomyocytes Recapitulating Calcium Handling Abnormalities. Front Physiol. 2018;9: 709. doi:10.3389/fphys.2018.00709
Tomek J, Bueno-Orovio A, Passini E, Zhou X, Minchole A, Britton O, et al. Development, calibration, and validation of a novel human ventricular myocyte model in health, disease, and drug block. eLife. 8: e48890. doi:10.7554/eLife.48890

## Repository structure 
**Primary folders and scripts**
- `GA`: contains main scripts to run genetic algorithm for model parameter calibration
- `Analysis`: contains scripts and functions to run post-calibration analyses
- `Pseudodataset`: contains scripts used to create the *in silico* Kernik model dataset used for protocol optimization
- `Classes_all`: code to construct classes using the corresponding cardiomyocyte model
- `dydts_all`: ODEs used to simulate electrophysiology dynamics based on each cardiomyocyte model
- `Other_helper_functions`: functions for setting up protocol simulations, reading in experimental or simulated data for the genetic algorithm, extracting and aligning action potential (AP) and calcium transient (CaT) traces, plotting calibrated parameter values, running independent validation simulations, etc.
- `Setup_and_parameters`: functions for setting up parameter values for some model classes in `Classes_all`

## Usage instructions
In MATLAB, open the project `CMmodelcalibration.prj` from the base folder.
**Generating a Kernik2019 model population and simulating APs & CaTs**
- *Note: to recreate the exact dataset used in the paper, skip steps 1-2 and start with step 3 using the conductance multipliers listed in Supplementary Table S2.*
1. Use `create_Kernik_cells` function in `Other_helper_functions` to initialize a population of *in silico* iPSC-CMs. Save `cells` to `Pseudodataset/saved_data/gaKernik_population_unfiltered.mat`. 
	- The non-default population parameters used in the paper were: `N = 100`, `sigma = 0.2`
2. Run `population_filtering.m` to filter the *in silico* population by spontaneous beating rate and randomly select a subset of model cells from the filtered population. 
3. Run `ga_pseudoDataset_protocols.m` to simulate AP and CaT traces for the chosen cells under the experimental conditions listed in the paper.
4. Run `ga_subset_processing.m` to process and save the *in silico* dataset to `Pseudodataset/saved_data`
   
**Running a single calibration to one dataset**
1. Before starting, make sure your AP and CaT datasets are in either of these folders: 
	- `Pseudodataset/saved_data` for *in silico* data
	- `ExperimentalData` for *in vitro* data
2. In the main script, `GA/sga_baseline_k19.m`: adjust the experiment information at the beginning of the script as needed, then run the entire script. Resulting calibrated parameter values and figures will be stored in `GA/Results`.
3. For scripts and instructions to run multiple calibration(s) on a job scheduler, please contact janice.yang@icahn.mssm.edu
   
**Analysis and evaluation of calibrated models**
- `Analysis/ga_fit_stats.m`: function to calculate calibration errors, calibration spreads, and other metrics from multiple GA runs
- `Analysis/aggregate_stats.m`: script to compile calibration error and spread metrics for multiple cells
- `plotFinalParams.m`: distribution plots of calibrated parameters from multiple GA runs
- `validation_threshold.m`: calculating predicted IKr block tolerance thresholds of calibrated models from GA runs
- Parameter sensitivity analysis scripts: to be added to repository
## Citation 
Yang J, Daily N, Pullinger TK, Wakatsuki T, Sobie EA. Creating cell-specific computational models of stem cell-derived cardiomyocytes using optical experiments. bioRxiv; 2024. p. 2024.01.07.574577. doi:10.1101/2024.01.07.574577

