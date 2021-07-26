function [results_struct, lo_table] = ...
    caramel_leaveOut(phenotype_structure, interaction_structure, varargin)

% DESCRIPTION: 
% This function conducts a leave-one-out analysis for training and
% testing CARAMeL models. 
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 2. Conduct leave-out analysis
% 3. Define output and display (if prompted)
% 
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
  
% I/O/U
%{
REQURED INPUTS: 
    1. phenotype_structure:     structure variable containing 3 fields
       a. data:         numeric matrix for drug profile data 
                        (e.g. chemogenomic, transcriptomic)
                        rows = features, columns = conditions
       b. features:     cell array of feature names (e.g. gene names)
       c. conditions:   cell array of conditions names (e.g. glucose)

    2. interaction_structure: 	structure variable containing 2 fields
                                'names' is required; choose between 
                                'scores' and 'class'
       a. names:        cell array of drug combination names, where
                        each row corresponds to a drug combination
          -> note: if list contains combinations of different lengths,
                   leave empty cells for combinations with smaller 
                   number of drugs
       b. scores:       numeric array of drug interaction scores that 
                        correspond to combinations in 'names'
       c. class:        cell or categorical array of interaction types 
                        that correspond to combinations in 'names'
          -> note: class names can only include 'Synergy', 'Additivity', 
                        and 'Antagonism' entries

OPTIONAL INPUTS: 
    1. LeaveOut:        string or char entry specifying what kind of
                        leave-out analysis to conduct
       -> default:      LeaveOut = 'all'
       -> note:         can also specify levels (e.g. 1, 2)

    2. ApplyLO:         boolean array specifying which entries in 
                        interaction data to apply LO analysis 
       -> default:      ApplyKO = true(size(scores))

    3. DisplayOutput:	boolean specifying whether to show results and 
                        plots for model performance
       -> default:      DisplayOutput = false)

OTHER FUNCTION PARAMETERS: 
    1. curdis:          any additional parameters available with the
                        'caramel' function
    
OUTPUTS: 
    1. lo_table:                table containing individual results for 
                                each unique drug in the interaction data
    2. predScores:              model predictions for all interactions

EXAMPLE USAGE: 
    1. Conduct a leave-all-out analysis given phenotype structure (PS) and 
       interaction structure (IS): 
       >> LO_table = caramel_leaveOut(PS, IS); 

    2. Conduct a leave-first-out analysis and return predictions: 
       >> [LO_table, predScores] = caramel_leaveOut(PS, IS, 'LeaveOut', 1); 

    3. Apply leave-second-out analysis for first 50 entries in IS and 
       display output table: 
       >> applylo = false(size(IS.scores)); 
       >> applylo(1:50) = true(50, 1);
       >> LO_table = caramel_leaveOut(PS, IS, 'ApplyLO', applylo, ...
                                     'DisplayOutput', true);
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define accepted ML functions and class types
    regML = {'regRF','fitrensemble','fitrsvm'};
    classML = {'classRF','fitcensemble','fitcsvm'}; 
    
    % Define required input validation functions
    isValidPhenotype = @(x) ...
        isstruct(x) && numel(fieldnames(x)) == 4 && ...
        all(isfield(x, {'data','features','conditions','H'})) && ...
        size(x.data, 1) == numel(x.features) && ...
        size(x.data, 2) == numel(x.conditions) && ...
        size(x.data, 2) == numel(x.H);
    isValidInteraction = @(x) ...
        isstruct(x) && isfield(x, 'names') && ...
        (isfield(x, 'scores') && size(x.names, 1) == numel(x.scores)) || ...
        (isfield(x, 'class') && size(x.names, 1) == numel(x.class));
    
    % Define required input variables
    addRequired(p, 'phenotype_structure', isValidPhenotype)
    addRequired(p, 'interaction_structure', isValidInteraction)
    
    % Define optional input validation functions 
    isValidLeaveOut = @(x) strcmpi(x, 'all') || ...
        (isnumeric(x) && numel(x) == 1 && rem(x, 1) == 0 && x > 0);
    isValidTimer = @(x) any(validatestring(x, {'on', 'off'}));
    isValidBoolean = @(x) islogical(x) || (isnumeric(x) && ...
        (x == 0 || x == 1));
    
    % Define optional input parameters
    addParameter(p, 'LeaveOut', 'all', isValidLeaveOut)
    addParameter(p, 'ApplyLO', [], isValidBoolean)
    addParameter(p, 'Timer', 'on', isValidTimer)
    addParameter(p, 'Verbose', false, isValidBoolean)
    addParameter(p, 'DisplayOutput', false, isValidBoolean)
    
    % Parse through inputs
    p.KeepUnmatched = true; 
    parse(p, phenotype_structure, interaction_structure, varargin{:})
    if ~isempty(fieldnames(p.Unmatched))
        extraParams = p.Unmatched; 
        if p.Results.Verbose
            disp('Extra name-value pair arguments provided:')
            disp(p.Unmatched)
        end
    else
        extraParams = [];
    end
    
    % Define list of optional input parameters from external functions
    caramelParams = {'Key', 'MLtype', 'MLmodel', 'Verbose', ...
        'regRFtrainOptions', 'classRFtrainOptions', 'classRFpredicOptions'};
    
    % Find any extra parameters and allocate to CARDIG function
    if ~isempty(extraParams)
        extraNames = fieldnames(extraParams); 
        [~, p_idx] = intersect(lower(extraNames), lower(caramelParams)); 
        if numel(p_idx) > 0
            caramelVarargin = cell(1, 2*numel(p_idx)); 
            for j = 1:numel(p_idx)
                [i1, i2] = deal(2*j-1, 2*j); 
                caramelVarargin{i1} = extraNames{p_idx(j)};
                eval(sprintf('caramelVarargin{i2} = extraParams.%s;', ...
                   extraNames{p_idx(j)}))
            end
            if strcmpi(caramelVarargin, 'MLtype')
                idx = find(strcmpi(caramelVarargin, 'MLtype')); 
                MLtype = caramelVarargin{idx + 1}; 
            else
                MLtype = 'regRF'; 
            end
        else
            caramelVarargin = {}; MLtype = 'regRF'; 
        end
    else
        caramelVarargin = {}; MLtype = 'regRF'; 
    end
    
    % Define inputs by name
    ps      = p.Results.phenotype_structure; 
    is      = p.Results.interaction_structure; 
    lo      = p.Results.LeaveOut; 
    applylo = p.Results.ApplyLO; 
    timer   = p.Results.Timer; 
    verbose = logical(p.Results.Verbose); 
    display = logical(p.Results.DisplayOutput); 
    
    % Start timer (if prompted)
    if strcmpi(timer, 'on')
        tStart = tic; 
    end
    
    % Display message
    if isnumeric(lo)
        s = num2str(lo); 
        switch s
            case '1'
                s = strcat(s, 'st'); 
            case '2'
                s = strcat(s, 'nd'); 
            case '3'
                s = strcat(s, 'rd'); 
            otherwise
                s = strcat(s, 'th'); 
        end
    elseif strcmpi(lo, 'all')
        s = lo; 
    else
        error('Invalid input provided for LeaveOut. Exiting function.')
    end
    if verbose
        fprintf('Initializing leave-%s-out analysis... \n', s)
    end

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check LeaveOut value
    if isnumeric(lo) && lo > size(is.names, 2)
        error('LeaveOut value is invalid.')
    end
    
    % Check ApplyLO value
    if isempty(applylo)
        applylo = true(size(is.names, 1), 1); 
    elseif size(applylo) ~= size(is.names, 1)
        error('Provide valid input for ApplyLO parameter.')
    end
    
    % Make sure required fields (based on MLtype) are present
    if ismember(MLtype, regML) && ~isfield(is, 'scores')
        error('Provide interaction scores for regression-based model.')
    end
    if ismember(MLtype, classML) && ~isfield(is, 'class')
        error('Provide interaction classes for classifier-based model.')
    end
    
    % Specify classification thresholds (if needed)
    if ~isfield(is, 'threshold')
        is.threshold = quantile(is.scores, [0.25 0.75]); 
    end
    
%% CONDUCT LEAVE-OUT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define unique list of conditions
    if isnumeric(lo)
        Condition = unique(is.names(applylo, lo)); 
    elseif strcmpi(lo, 'all')
        Condition = unique(is.names(applylo, :)); 
    end
    Condition(cellfun(@isempty, Condition)) = []; 
    n = numel(Condition); 
    
    % Leave-out analysis
    predScores = nan(sum(applylo), 1); predClass = cell(sum(applylo), 1); 
    [N, R, p, Accuracy] = deal(nan(n, 1), nan(n, 1), nan(n, 1), nan(n, 1)); 
    [synAUC, antAUC] = deal(nan(n, 1), nan(n, 1)); 
    progressbar(sprintf('Conducting leave-%s-out analysis...', s))
    for i = 1:n
        % Determine train and test indices
        if isnumeric(lo)
            idx = applylo & strcmpi(is.names(:, lo), Condition{i}); 
        elseif strcmpi(lo, 'all')
            idx = applylo & sum(ismember(is.names, Condition{i}), 2) > 0; 
        end
        % Define training structure
        train_struct = struct(...
            'names', {is.names(~idx, :)}, 'scores', is.scores(~idx)); 
        % Define testing structure
        test_struct = struct(...
            'names', {is.names(idx, :)}, 'scores', is.scores(idx), ...
            'threshold', is.threshold); 
        % Train caramel model
        [~, model] = caramel(ps, train_struct, 'train', caramelVarargin{:}); 
        % Predict interaction scores using GEM-ML model
        data = caramel(ps, test_struct, 'predict', ...
            'MLmodel', model, caramelVarargin{:});
        results = caramel_assessPerf(data, 'InteractionType', 'Synergy'); 
        % Save outputs of interest
        ix = find(idx(applylo)); 
        predScores(ix) = results.predScores; 
        predClass(ix) = results.predClass; 
        N(i) = size(results.interactionNames, 1); 
        R(i) = results.corrStat; 
        p(i) = results.corrPvalue; 
        Accuracy(i) = results.accuracy; 
        synAUC(i) = results.synergyAUC; 
        antAUC(i) = results.antagonismAUC; 
        % Update progress bar
        progressbar(i/n)
    end
    
%% DEFINE OUTPUT AND DISPLAY (IF PROMPTED) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define outputs
    data_struct = caramel_processInteractions(is, 'Indices', applylo); 
    data_struct.predScores = predScores; 
    data_struct.predClass = predClass;
    data_struct.MLmethod = data.MLmethod; 
    results_struct = caramel_assessPerf(data_struct); 
    lo_table = table(Condition, N, R, p, Accuracy, synAUC, antAUC); 
    
    % Display output (if prompted)
    if display
        disp(lo_table)
    end
    if verbose
        fprintf('Leave-%s-out analysis completed. \n', s)
    end
    
    % End timer (if needed)
    if strcmpi(timer, 'on')
        tEnd = toc(tStart);
        [m, s] = deal(floor(tEnd/60), round(rem(tEnd, 60)));
        fprintf('Elapsed time is %d minute(s) and %d second(s) \n', m, s)
    end
    
end