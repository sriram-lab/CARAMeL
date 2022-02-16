function flux_struct = derive_flux(model, gene_struct, varargin)

    % DESCRIPTION 
    % This function derives simulated steady-state fluxes through metabolic
    % reactions given a genome-scale metabolic model (GEM) and differential
    % gene expression data (or similar information). 
    % 
    % STEPS 
    % 1. Parse through inputs
    % 2. Ascertain input compatibility
    % 3. Determine differentially expressed genes (DEGs)
    % 4. Derive simulated reaction fluxes
    % 5. Return flux data and save (if prompted)
    % 
    % Author:   Carolina H. Chung
    % Created:  September 30, 2019
    % Updates:  September 20, 2020 (Carolina H. Chung)
    
    % I/O
    %{
    REQUIRED INPUTS: 
        1. model:       genome-scale metabolic model
           -> type:     MATLAB structure
    
        2. gene_struct: information on differentially expressed genes
           -> type:     MATLAB structure
           -> dim:      m x 1 (m = number of fields; at least 3)
           -> note:     must contain at least the following 3 fields
                        1. genes:       cell array of gene names
                        2. conditions:  cell array of treatment groups
                        3. fold_data:   double array of gene expression 
                                        fold change
                        -> length of 'genes' and 'conditions' must match
                           row and column length in 'fold_data'
           -> note:     may also contain a field called 'sig_data'
                        if not specified, all fold values are assumed to be
                        significant
    
    OPTIONAL INPUTS:
        1. ModelSpecs:  additional constraints to use for simulating flux
                        (i.e. kappa, rho, epsilon)
           -> type:     integer or double vector
           -> dim:      1 x 3 (where 1 = kappa, 2 = rho, 3 = epsilon)
           -> note:     refer to `constrain_flux_regulation.m` for more
                        details on each parameter 
    
        2. Threshold:   threshold values defining down- and up-regulation
           -> type:     integer or double vector 
           -> dim:      1 x 2 (1st element = down threshold, 
                               2nd element = up threshold)
    
        3. Reverse:     Boolean indicating if threshold convention should
                        be reversed
           -> note:     set to true if using chemogenomic data
    
        4. SaveData:    Boolean indicating whether output should be saved
                        into an Excel file
    
        5. Filename:    file to save simulated flux data into
           -> type:     string/char
    
        6. WRITETABLE:  any additional inputs for 'writetable' function can
                        be included as Name-Value pairs when writing to
                        'filename'
    
    OUTPUTS:
        1. flux_struct: MATLAB structure containing the following fields 
           -> degs:     MATLAB table of DEGs for each condition
           -> grRate:   simulated growth rate for each condition
           -> flux:     MATLAB table of simulated flux for each condition;
                        first 3 columns list reaction IDs, reactions names,
                        and sub-systems
    
    EXAMPLE USAGE: 
        1. Default usage: 
           >> flux_struct = derive_flux(model, gene_struct)
    
        2. Specify threshold and save into 'flux_output.xlsx' with 
           sheet name 'ecoli':  
           >> flux_struct = derive_flux(model, gene_struct, ...
                                        'Threshold', [-2 2], ...
                                        'SaveData', true, ...
                                        'Sheetname, 'ecoli')
    %}
    
%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidGeneStruct = @(x) isstruct(x) && ...
        all(ismember({'genes','conditions','fold_data'}, fieldnames(x))); 
    isValidModelSpecs = @(x) isnumeric(x) && (numel(x) == 3); 
    isValidThreshold = @(x) isnumeric(x) && (numel(x) == 2); 
    isValidBoolean = @(x) islogical(x) || ismember(x, [0 1]); 
    isValidName = @(x) isstring(x) || ischar(x); 
    
    % Define required input variables
    addRequired(p, 'model', @isstruct)
    addRequired(p, 'gene_struct', isValidGeneStruct)
    
    % Define optional input parameters
    addParameter(p, 'ModelSpecs', [1E-2 1E-2 1], isValidModelSpecs)
    addParameter(p, 'Threshold', [-2 2], isValidThreshold)
    addParameter(p, 'Reverse', false, isValidBoolean)
    addParameter(p, 'SaveData', false, isValidBoolean)
    addParameter(p, 'Filename', 'flux_data.xlsx', isValidName)
    addParameter(p, 'Verbose', false, @islogical)
    
    % Parse through inputs
    p.KeepUnmatched = true; 
    parse(p, model, gene_struct, varargin{:})
    
    % Define inputs by name
    model           = p.Results.model; 
    gene_struct     = p.Results.gene_struct; 
    modelSpecs      = p.Results.ModelSpecs; 
    threshold       = p.Results.Threshold; 
    reverse         = p.Results.Reverse; 
    saveData        = p.Results.SaveData; 
    filename        = p.Results.Filename; 
    verbose         = p.Results.Verbose;
    
    % Define list of optional input parameters from external functions
    writetableParams = {'FileType','WriteVariableNames',...
        'WriteRowNames','DateLocale','WriteMode','Delimiter',...
        'QuoteStrings','Encoding','Sheet','Range','UseExcel'}; 
    
    % Find any optional parameters and allocate to appropriate function
    if ~isempty(fieldnames(p.Unmatched))
        disp('Extra name-value pair arguments provided:')
        disp(p.Unmatched)
        extraParams = p.Unmatched; 
        extraNames = fieldnames(extraParams); 
        % WRITETABLE function
        [~, p_idx] = intersect(extraNames, writetableParams);
        if numel(p_idx) > 0
            writetableVarargin = cell(1, 2*numel(p_idx));
            writetableVarargin(1:2:end-1) = extraNames(p_idx);
            i = 0; 
            for j = 2:2:length(writetableVarargin)
                i = i + 1; 
                writetableVarargin(j) = ...
                    {eval(sprintf('extraParams.%1$s', extraNames{i}))};
            end
            disp(['Extra name-value pair argument(s) passed to ', ...
                'WRITETABLE function.'])
        else
            writetableVarargin = {};
        end
    end

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check model
    model_fields = {'genes','geneNames','rxns','rxnNames','subSystems'};
    if ~all(ismember(model_fields, fieldnames(model)))
        error(['At least one of the following fields is not present ' ...
            'in model: ' strjoin(model_fields)])
    end
    
    % Check gene structure
    if numel(gene_struct.genes) ~= size(gene_struct.fold_data, 1)
        error("Dimension mismatch between 'genes' and 'fold_data'.")
    end
    if numel(gene_struct.conditions) ~= size(gene_struct.fold_data, 2)
        error("Dimension mismatch between 'conditions' and 'fold_data'.")
    end
    if ~ismember('sig_data', fieldnames(gene_struct))
        gene_struct.sig_data = 0.01 .* ones(size(gene_struct.fold_data)); 
    end
    
    % Check that genes in gene_struct match either model genes or geneNames
    if ~any(ismember(gene_struct.genes, model.genes)) && ...
            ~any(ismember(gene_struct.genes, model.geneNames))
        error(strcat("Genes in 'gene_struct' do not match either ", ...
            "model genes or geneNames."))
    end
    
    % Start function timer
    if verbose
        tStart = tic; 
    end

%% DETERMINE DIFFERENTIALLY EXPRESSED GENES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create empty variables to store DEGs
    Down = cell(numel(gene_struct.conditions), 1); 
    Up = cell(size(Down)); 
    
    % Iterate through each condition (i.e. treatment group)
    for j = 1:length(gene_struct.conditions)
%         threshold = quantile(gene_struct.fold_data(:, j), [0.2 0.8]); 
        if reverse
            up_idx = (gene_struct.sig_data(:, j) < 0.05) & ...
                (gene_struct.fold_data(:, j) < threshold(1));
            down_idx = (gene_struct.sig_data(:, j) < 0.05) & ...
                (gene_struct.fold_data(:, j) > threshold(2));
        else
            down_idx = (gene_struct.sig_data(:, j) < 0.05) & ...
                (gene_struct.fold_data(:, j) < threshold(1));
            up_idx = (gene_struct.sig_data(:, j) < 0.05) & ...
                (gene_struct.fold_data(:, j) > threshold(2));
        end
        if any(ismember(gene_struct.genes, model.geneNames))
            [~, ~, down_mod_idx] = ...
                intersect(gene_struct.genes(down_idx), model.geneNames);
            [~, ~, up_mod_idx] = ...
                intersect(gene_struct.genes(up_idx), model.geneNames); 
        else
            [~, ~, down_mod_idx] = ...
                intersect(gene_struct.genes(down_idx), model.genes);
            [~, ~, up_mod_idx] = ...
                intersect(gene_struct.genes(up_idx), model.genes); 
        end
        Down{j} = model.genes(down_mod_idx); 
        Up{j} = model.genes(up_mod_idx); 
    end
    if verbose
        disp('Up- and down-regulated genes determined for all conditions.')
    end

%% DERIVE SIMULATED REACTION FLUXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define empty variables to store output data
    flux = zeros(length(model.rxns), length(gene_struct.conditions)); 
    growth = zeros(length(gene_struct.conditions), 1); 
    
    % Simulate flux activity based on DEGs
    count = 0;
    kappa = modelSpecs(1); rho = modelSpecs(2); epsilon = modelSpecs(3); 
    if verbose
        progressbar('Deriving reaction fluxes...')
    end
    for j = 1:length(gene_struct.conditions) 
        try
            [flux(:, j), growth(j), ~] = ...
                constrain_flux_regulation(model, Up{j}, Down{j}, ...
                kappa, rho, epsilon, 0);
        catch
            obj_pos = model.c; 
            [flux(:, j), ~, ~] = ...
                constrain_flux_regulation(model, Up{j}, Down{j}, ...
                kappa, rho, epsilon, 0);
            growth(j) = flux(obj_pos, j); 
        end
        if isempty(Up{j}) && isempty(Down{j})
            warning(['No DEGs found for ' gene_struct.conditions{j} ...
                '. Inputting default flux and growth values.'])
            count = count + 1; 
        end
        if verbose
            progressbar(j / length(gene_struct.conditions))
%             if reverse
%                 progressbar([], j / length(gene_struct.conditions))
%             else
%                 progressbar([], [], j / length(gene_struct.conditions))
%             end
        end
    end
    if verbose
        disp('Growth and flux determined.')
        fprintf('Total warnings: %d out of %d\n', count, j);
    end

%% RETURN FLUX DATA AND SAVE (IF PROMPTED) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create flux table
    rxn_table = table(model.rxns, model.rxnNames, model.subSystems, ...
        'VariableNames', {'rxn','rxnName','subSystem'}); 
    flux_table = array2table(flux, 'VariableNames', gene_struct.conditions');
    final_table = horzcat(rxn_table, flux_table); 
    
    % Define output 
    flux_struct = struct('degs', {table(Up, Down)}, ...
        'grRate', growth, 'flux', {final_table}); 
    
    % Save data (if prompted)
    if saveData
        writetable(which(filename), writetableVarargin{:})
        disp(['Location: ' which(filename)])
    end
    if verbose
        tEnd = toc(tStart);
        fprintf('Elapsed time is %d minutes and %f seconds.\n', ...
            floor(tEnd / 60), rem(tEnd, 60));
    end
    
end