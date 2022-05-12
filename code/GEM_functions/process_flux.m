function phenotype_struct = process_flux(flux_table, model, varargin)

    % DESCRIPTION 
    % This function process flux data into Boolean-based phenotype data
    % indicative of differential flux activity. 
    % 
    % STEPS 
    % 1. Parse through inputs
    % 2. Ascertain input compatibility
    % 3. Process data in flux table
    % 4. Determine differential flux activity
    % 5. Define output
    % 
    % Author:   Carolina H. Chung
    % Created:  September 21, 2020
    % Updates:  September 21, 2020 (Carolina H. Chung)
    
    % I/O
    %{ 
    REQUIRED INPUTS: 
        1. flux_table:  MATLAB table of GEM-based flux data
           -> type:     MATLAB table
           -> dim:      m x n + 3 (m = # of reactions in GEM, 
                                   n = # of treatment groups (conditions))
           -> note:     first 3 columns list reactions, reaction names, and
                        corresponding sub-systems
           -> note:     flux output from `derive_flux.m`
    
        2. Model:       GEM from which flux data was simulated
           -> note:     provide if 'Control' is not specified
           -> default:  model = struct()   

    OPTIONAL INPUTS:     
        1. Threshold:   threshold values defining differential flux
           -> type:     integer or double vector 
           -> dim:      1 x 2 (1st element = negative threshold, 
                               2nd element = positive threshold)
           -> default:  Threshold = [-2 2]
    
        2. Label:       specify what info to use for 'phenotype_label'
           -> type:     string or char
           -> note:     choose from {'rxn','rxnName','subSystem'}
           -> default:  Label = 'rxn'
    
        3. RemoveZeros: Boolean for removing all-zero vectors
           -> type:     Boolean
           -> default:  RemoveZeros = true

    OUTPUTS: 
        1. phenotype_struct:    MATLAB structure containing the following
           a. data:             Boolean matrix, where entries of 1 are 
                                associated with differential flux
              -> type:          numerical matrix of 0's and 1's
              -> dim:           m x n (m = length of phenotype_labels, 
                                       n = length of conditions)
           b. features:         reaction names corresponding to data rows
              -> type:          cell array
              -> dim:           m x 1
    
           c. conditions:       treatment group labels corresponding to
                                data columns
              -> type:          cell array
              -> dim:           1 x n
    
    EXAMPLE USAGE: 
        1. Default usage: 
           >> flux_struct = process_flux(flux_table, model)
    
        2. Specify threshold = 1.5: 
           >> flux_struct = process_flux(flux_table, model, ...
                                         'Threshold, 1.5)
    %}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidFluxTable = @(x) istable(x) && ...
        all(ismember({'rxn','rxnName','subSystem'}, fieldnames(x))); 
    isValidThreshold = @(x) isnumeric(x) && (numel(x) == 2) && ...
        x(1) < 0 && x(2) > 0; 
    isValidName = @(x) isstring(x) || ischar(x) && ...
        ismember(x, {'rxn','rxnName','subSystem'}); 
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 0 || x == 1)); 
    
    % Define required input variables
    addRequired(p, 'flux_table', isValidFluxTable)
    addRequired(p, 'model', @isstruct)
    
    % Define optional input parameters
    addParameter(p, 'Threshold', [-2 2], isValidThreshold)
    addParameter(p, 'Label', 'rxn', isValidName)
    addParameter(p, 'RemoveZeros', false, isValidBoolean)
    
    % Parse through inputs
    parse(p, flux_table, model, varargin{:})
    
    % Define inputs by name
    flux_table  = p.Results.flux_table; 
    model       = p.Results.model; 
    threshold   = p.Results.Threshold; 
    label       = p.Results.Label; 
    remZs       = logical(p.Results.RemoveZeros); 

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check that flux_table metadata matches model fields
    if eval(sprintf('~all(strcmp(flux_table.%1$s, model.rxns))', label))
        error(['Specified feature labels do not match corresponding ' ...
            'data in flux table.'])
    end
%     if ~all(strcmp(flux_table.rxn, model.rxns))
%         error('Flux table reactions do not match model reactions.')
%     end
%     if ~all(strcmp(flux_table.rxnName, model.rxnNames))
%         error(['Flux table reaction names do not match model ' ...
%             'reaction names.'])
%     end
%     if ~all(strcmp(flux_table.subSystem, model.subSystems))
%         error(['Flux table reaction sub-systems do not match model ' ...
%             'reaction sub-systems.'])
%     end
    
%% PROCESS DATA IN FLUX TABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Extract data of interest
    flux = table2array(flux_table(:, 4:end)); 
    if iscell(flux)
        flux = cell2mat(flux); 
    end
    labels = eval(sprintf('flux_table.%1$s', label)); 
    
    % Simulate flux for control (WT)
    z = constrain_flux_regulation(model, {}, {}, [], [], [], 0); 
    
    % Remove all-zero rows (if prompted)
    if remZs
        zero_idx = sum([flux z], 2) == 0; 
        flux(zero_idx,:) = []; labels(zero_idx) = []; z(zero_idx) = [];
    end
    
    % Convert zeros to small non-zero value
    abs_num = abs([flux(:); z(:)]); 
    abs_min = min(abs_num(abs_num ~= 0), [], 'all'); 
    flux(~flux) = 1E-3 * abs_min; z(~z) = 1E-3 * abs_min;
    
    % Normalize data
    data = (flux - z) ./ z; 

%% DETERMINE DIFFERENTIAL FLUX ACTIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Negative flux change (compared to WT)
    neg_flux = data < threshold(1);
    neg_labels = strcat('neg-', labels); 
    if remZs
        neg_idx = sum(neg_flux, 2) == 0; 
        neg_flux(neg_idx, :) = []; 
        neg_labels(neg_idx) = [];
    end
    
	% Positive flux change (compared to WT)
    pos_flux = data > threshold(2); 
    pos_labels = strcat('pos-', labels);
    if remZs
        pos_idx = sum(pos_flux, 2) == 0; 
        pos_flux(pos_idx, :) = []; 
        pos_labels(pos_idx) = [];
    end

    % Metabolic flux entropy
    H = log(var(flux)); 
    
%% DEFINE OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    phenotype_struct = struct('data', [neg_flux; pos_flux], 'H', H, ...
        'features', {vertcat(neg_labels, pos_labels)}, ...
        'conditions', {flux_table.Properties.VariableNames(4:end)}); 

end