# CARAMeL
Condition-specific Antibiotic Regimen Assessment using Mechanistic Learning

# License
Released via GPL GNU License  
&copy; 2022 The Regents of the University of Michigan  
Chandrasekaran Research Group - https://systemsbiologylab.org/  
Contact: csriram@umich.edu  

# Summary
This GitHub repository contains all code files used in data analysis and model development for the CARAMeL publication (under review). 

# Dependencies
Successful implementation of the code files within this repository require the following external installations: 
- [MATLAB](https://www.mathworks.com/products/matlab.html): preferrably version 2019b or higher
- [Gurobi](https://www.gurobi.com/): optimization solver; preferrably latest version
- [COBRA Toolbox](https://github.com/opencobra/cobratoolbox): preferrably clone directly from GitHub

# Installation and execution: 
1. Install a local version of this repository using the following command in a terminal: 
```
git clone https://github.com/sriram-lab/CARAMeL.git
```
2. Ensure to add the CARAMeL repository into your MATLAB path (see instructions [here](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html)). 
3. You should now be able to run all MATLAB livescript files in the main repository folder.

# Repository structure
The repository structure is outlined below. Of note, **folders** are designated in bold while *data/code files* are designated in italics. 
```
CARAMeL
└───**code**: contains all relevant code files
└──────**CARAMeL_suite**: contains all code files created to implement the CARAMeL approach. Note: all files contain details on inputs, outputs, and function usage
│   │   │   *caramel.m*:                      main script the constructs a CARAMeL model
│   │   │   *caramel_anova.m*:                conducts a one-way ANOVA test to determine GEM reactions with differential flux activity 
│   │   │   *caramel_assessPerf.m*:           assesses the predictive performance of a CARAMeL model
│   │   │   *caramel_auroc*:                  calculates the area under the receiver operating curve (AUROC) for a CARAMeL model 
│   │   │   *caramel_classify*:               classifies drug interactions as synergistic, additive, or antagonistic given interaction scores and threshold information
│   │   │   *caramel_crossValidation.m*:      conducts a k-fold cross-validation for a CARAMeL model
│   │   │   *caramel_features2gem.m*:         extracts GEM reactions associated with top CARAMeL model features
│   │   │   *caramel_leaveOut.m*:             conducts a leave-out analysis for a CARAMeL model
│   │   │   *caramel_plot.m*:                 generates plots visualizing CARAMeL model performance
│   │   │   *caramel_processInteractions.m*:  processes drug interaction data into a standardized format
│   │   │   *caramel_rankSubsystems.m*:       determines significantly enriched metabolic pathways based on GEM reactions associated with top CARAMeL model features
│   │   │   *caramel_screen.m*:               screens all possible drug combination outcomes given a list of drugs and combination order
│   │   │   *caramel_topFeatures.m*:          extracts the top CARAMeL model features 
└──────**GEM_functions**: contains all code files relevant to simulating metabolism using GEMs
│   │   │   *change_media.m*:                 changes the simulated media condition for a given GEM
│   │   │   *constrain_flux_regulation.m*:    simulates reaction fluxes based on omics data constraints
│   │   │   *derive_flux.m*:                  derives flux simulations for a specified list of conditions
│   │   │   *process_flux.m*:                 processes simulated flux data into phenotypic data to be used for CARAMeL model construction
└──────**misc**: folder containing miscellaneous code files
└─────────**File_Exchange**: code available through MATLAB File Exchange. Direct download links are provided for each method
└────────────**confusion matrix**:            determines the confusion matrix and other relevant information for two or more classes ([download](https://www.mathworks.com/matlabcentral/fileexchange/60900-multi-class-confusion-matrix))
└────────────**hline_vline**:                 generates horizontal and vertical lines on MATLAB plots ([download](https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline))
└────────────**multiple_boxplot**:            generates multiple boxplots in a single MATLAB figure ([download](https://www.mathworks.com/matlabcentral/fileexchange/47233-multiple_boxplot-m))
└────────────**permn.m**:                     determines all possible permutations with repetition ([download](https://www.mathworks.com/matlabcentral/fileexchange/7147-permn))
└────────────**progressbar**:                 displays a progress bar to run a piece of code ([download](https://www.mathworks.com/matlabcentral/fileexchange/6922-progressbar))
└────────────**raacampbell-sigstar-5aabaeb**: plots statistical significance stars onto MATLAB plots ([download](https://www.mathworks.com/matlabcentral/fileexchange/39696-raacampbell-sigstar))
└─────────**custom**: custom code files created for specific purposes. Refer to function file for details on input/output/usage
│   │   │   │   *extract_from_confMatStats.m*:    extracts information from a confusion matrix
│   │   │   │   *extract_metadata.m*:             extracts metadata from an Excel file
│   │   │   │   *gscatter2.m*:                    generates a 2D group scatter plot
│   │   │   │   *meansgraph.m*:                   customized modification to the built-in MATLAB function of the same name
│   │   │   │   *process_chemgen_v2.m*:           processes chemogenomic data for *E. coli* ([source](https://doi.org/10.15252/msb.20156777))
│   │   │   │   *process_transcriptome_tb.m*:     processes transcriptomic data for *M. tb* ([source](https://doi.org/10.1128/mBio.02627-19))
│   │   │   │   *sig_boxplot.m*:                  generates a multple boxplot figure with significance annotations
└───**data**: contains all data files used for CARAMeL model development and downstream analyses
└──────**supplementary**: contains all supplementary data files
│   │   │   *ecoli_proteomics.xlsx*:          proteomics data for *E. coli*
│   │   │   *ecoli_transcriptomics.xlsx*:     transcriptomics data for *E. coli*
│   │   *biolog_PM01_conditions.xlsx*:        list of all conditions in the Biolog phenotype microarray 1 (PM01) ([source](https://www.biolog.com/products-portfolio-overview/phenotype-microarrays-for-microbial-cells/))
│   │   *ecoli_chemogenomics.xlsx*:           chemogenomic data for *E. coli* ([source](https://doi.org/10.1016/j.cell.2010.11.052))
│   │   *ecoli_interactions.xlsx*:            drug interaction data for *E. coli*
│   │   *ecoli_media.xlsx*:                   details on different media conditions to simulate for the *E. coli* GEM
│   │   *iEK1008.mat*:                        *M. tb* GEM ([source](https://doi.org/10.1186/s12918-018-0557-y))
│   │   *iJO1366.xml*:                        *E. coli* GEM ([source](https://dx.doi.org/10.1038%2Fmsb.2011.65))
│   │   *mtb_interactions.xlsx*:              drug interaction data for *M. tb*
│   │   *mtb_transcriptomics.xlsx*:           transcriptomic data for *M. tb* ([source](https://doi.org/10.1128/mBio.02627-19))
└───**results**: contains all results from CARAMeL model analyses
│   │   *CARAMeL_summary.xlsx*:               summary data from the *E. coli* CARAMeL model
│   │   *plot_data.xlsx*:                     Excel file containing data pertinent for figure generation
│   │   *CARAMeL_workspace.mat*:              [REMOVED; CONTACT chechung@umich.edu FOR FILE] MATLAB workspace file containing results for CARAMeL
│   *CARAMeL_figures.mlx*:        MATLAB livescript file that generates plots used in CARAMeL figures
│   *CARAMeL_flux.mlx*:           MATLAB livescript file that generates GEM flus simulation data for *E. coli* and *M. tb* CARAMeL models
│   *CARAMeL_main.mlx*:           MATLAB livescript file that conducts all main analyses for CARAMeL approach
│   *CARAMeL_validations.mlx*:    MATLAB livescript file that executes validations discussed in CARAMeL supplementary materials
│   *README.md*:                  markdown file for CARAMeL readme documentation
│   *data_visualization.Rmd*:     R markdown file generating additional plots used for CARAMeL figures
│   *data_visualization.html*:    HTML file of the R markdown file
```
