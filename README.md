# CARAMeL
Condition-specific Antibiotic Regimen Assessment using Mechanistic Learning

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
- **code**: 
  - **CARAMeL_suite**: 
    - *caramel.m*: 
    - *caramel_anova.m*:
    - *caramel_assessPerf.m*: 
    - *caramel_auroc*: 
    - *caramel_classify*: 
    - *caramel_crossValidation.m*: 
    - *caramel_features2gem.m*: 
    - *caramel_leaveOut.m*: 
    - *caramel_plot.m*: 
    - *caramel_processInteractions.m*: 
    - *caramel_rankSubsystems.m*: 
    - *caramel_screen.m*: 
    - *caramel_topFeatures.m*: 
  - **GEM_functions**: 
    - *change_media.m*: 
    - *constrain_flux_regulation.m*: 
    - *derive_flux.m*: 
    - *process_flux.m*: 
  - **misc**: 
    - **File_Exchange**: 
      - **confusion matrix**: 
      - **hline_vline**: 
      - **multiple_boxplot**: 
      - **permn.m**: 
      - **progressbar**: 
      - **raacampbell-sigstar-5aabaeb**: 
      - **randomforest-matlab**: 
    - **custom**: 
      - *extract_from_confMatStats.m*: 
      - *extract_metadata.m*: 
      - *gscatter2.m*: 
      - *meansgraph.m*: 
      - *process_chemgen_v2.m*: 
      - *process_transcriptome_tb.m*: 
      - *sig_boxplot.m*: 
- **data**: 
  - **supplementary**: 
    - *ecoli_proteomics.xlsx*:
    - *ecoli_transcriptomics.xlsx*: 
  - *biolog_PM01_conditions.xlsx*: 
  - *ecoli_chemogenomics.xlsx*: 
  - *ecoli_interactions.xlsx*: 
  - *ecoli_media.xlsx*: 
  - *iEK1008.mat*: 
  - *iJO1366.xml*: 
  - *mtb_interactions.xlsx*: 
  - *mtb_transcriptomics.xlsx*: 
- **results**: 
  - *CARAMeL_summary.xlsx*: 
  - *CARAMeL_workspace.mat*: 
  - *plot_data.xlsx*: 
- *CARAMeL_figures.mlx*: 
- *CARAMeL_flux.mlx*: 
- *CARAMeL_main.mlx*:
- *CARAMeL_validations.mlx*:
- *data_visualization.Rmd*: 
- *data_visualization.html*: 
